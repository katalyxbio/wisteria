use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::mpsc::sync_channel;
use std::time::{Duration, Instant};
use crate::debug::IoStats;
use crate::stats::QCAccumulator;
use needletail::{parse_fastx_file, parse_fastx_reader};
use noodles::bgzf;
use rayon::prelude::*;
use std::io::{Read, Write};

type FastqBatch = Vec<(Vec<u8>, Vec<u8>)>;

const BATCH_SIZE: usize = 10_000;
/// Batches the reader may run ahead of the folder; bounds memory and provides
/// backpressure so we never buffer the whole file.
const CHANNEL_DEPTH: usize = 2;

/// Detect whether `path` is BGZF (block-gzip, as written by `bgzip`/htslib)
/// rather than a plain single-stream gzip.
///
/// BGZF is still valid gzip, but it's a concatenation of independent ≤64 KB
/// blocks carrying a `BC` extra subfield — which is exactly what lets us
/// decompress it in parallel (see [`process_fastq`]). A plain `.fastq.gz` is
/// one continuous DEFLATE stream and can only be decoded serially.
///
/// We check the gzip magic + DEFLATE method + FEXTRA flag, then the `BC`
/// subfield id at the fixed offset `bgzip` writes it. Returns `Ok(false)` for
/// short files, plain gzip, or uncompressed input.
fn is_bgzf(path: &Path) -> std::io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 18];
    if let Err(e) = file.read_exact(&mut buf) {
        // Too short to be BGZF (which always has an 18-byte header) → not BGZF.
        if e.kind() == std::io::ErrorKind::UnexpectedEof {
            return Ok(false);
        }
        return Err(e);
    }
    Ok(buf[0] == 0x1f            // gzip magic
        && buf[1] == 0x8b
        && buf[2] == 0x08        // CM = DEFLATE
        && (buf[3] & 0x04) != 0  // FLG.FEXTRA
        && buf[12] == b'B'       // BGZF subfield id: 'B'
        && buf[13] == b'C')      // 'C'
}

pub fn process_fastq(
    path: &Path,
    decompress_threads: usize,
    min_len: Option<usize>,
    min_qual: Option<f32>,
    debug_stats: bool,
) -> Result<(QCAccumulator, IoStats), Box<dyn std::error::Error>> {
    // A plain `.fastq.gz` is one serial DEFLATE stream — the sole bottleneck on
    // large files. If the input is instead BGZF (block-gzip from `bgzip`), route
    // it through the same multithreaded block decompressor as the BAM path so
    // decompression parallelizes across cores; needletail then parses the
    // already-decompressed plaintext (it adds no gzip layer when it sees no
    // magic). Otherwise fall back to needletail's own (serial) gzip handling.
    let bgzf_parallel = is_bgzf(path).unwrap_or(false);

    // needletail's `Box<dyn FastxReader>` is not `Send`, so it stays on this
    // (producer) thread; the consumer thread owns the accumulator. Only the
    // owned `(seq, qual)` byte batches cross the channel, overlapping the
    // decompress/parse with the statistics fold.
    let mut reader = if bgzf_parallel {
        let file = File::open(path)?;
        let workers =
            NonZeroUsize::new(decompress_threads).unwrap_or(NonZeroUsize::new(1).unwrap());
        let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(workers, file);
        parse_fastx_reader(bgzf_reader)?
    } else {
        parse_fastx_file(path)?
    };

    let (tx, rx) = sync_channel::<FastqBatch>(CHANNEL_DEPTH);

    let consumer = std::thread::spawn(move || {
        let mut accumulator = QCAccumulator::new();
        let mut proc_time = Duration::default();
        let mut merge_time = Duration::default();
        let mut recv_wait = Duration::default();
        let mut batches = 0u64;
        let mut total_processed = 0usize;

        loop {
            // Time spent here is the folder sitting idle waiting for the reader.
            let recv_start = Instant::now();
            let batch = match rx.recv() {
                Ok(batch) => batch,
                Err(_) => break,
            };
            recv_wait += recv_start.elapsed();

            let proc_start = Instant::now();
            let local_acc = process_fastq_batch(&batch, min_len, min_qual);
            proc_time += proc_start.elapsed();

            let merge_start = Instant::now();
            total_processed += batch.len();
            accumulator.merge(local_acc);
            merge_time += merge_start.elapsed();
            batches += 1;

            if !debug_stats {
                print!("\rProcessed {} reads...", total_processed);
                let _ = std::io::stdout().flush();
            }
        }

        if !debug_stats && total_processed > 0 {
            println!("\rProcessed {} reads. Done!", total_processed);
        }

        (accumulator, proc_time, merge_time, recv_wait, batches)
    });

    let mut total_io_time = Duration::default();
    let mut send_wait = Duration::default();
    let mut read_error: Option<Box<dyn std::error::Error>> = None;
    'outer: loop {
        let io_start = Instant::now();
        let mut batch: FastqBatch = Vec::with_capacity(BATCH_SIZE);
        let mut ended = false;
        for _ in 0..BATCH_SIZE {
            match reader.next() {
                Some(Ok(seqrec)) => {
                    let seq = seqrec.seq().to_vec();
                    let qual = seqrec.qual().map(|q| q.to_vec()).unwrap_or_default();
                    batch.push((seq, qual));
                }
                Some(Err(e)) => {
                    read_error = Some(Box::new(e));
                    ended = true;
                    break;
                }
                None => {
                    ended = true;
                    break;
                }
            }
        }
        total_io_time += io_start.elapsed();

        if !batch.is_empty() {
            // Time spent here is the reader blocked because the folder is behind.
            let send_start = Instant::now();
            let sent = tx.send(batch);
            send_wait += send_start.elapsed();
            if sent.is_err() {
                break 'outer;
            }
        }
        if ended {
            break;
        }
    }
    drop(tx); // close the channel so the consumer's loop terminates

    let (accumulator, proc_time, merge_time, recv_wait, batches) = consumer
        .join()
        .map_err(|_| Box::<dyn std::error::Error>::from("FASTQ consumer thread panicked"))?;

    if let Some(e) = read_error {
        return Err(e);
    }

    let stats = IoStats {
        io_time: total_io_time,
        proc_time,
        merge_time,
        batches,
        send_wait,
        recv_wait,
        bgzf_parallel,
    };

    Ok((accumulator, stats))
}

fn process_fastq_batch(
    batch: &[(Vec<u8>, Vec<u8>)],
    min_len: Option<usize>,
    min_qual: Option<f32>,
) -> QCAccumulator {
    batch
        .par_iter()
        .fold(QCAccumulator::new, |mut acc, (seq, qual)| {
            acc.add_read(seq, qual, true, min_len, min_qual);
            acc
        })
        .reduce(QCAccumulator::new, |mut a, b| {
            a.merge(b);
            a
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTQ: &[u8] = b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACG\n+\nIII\n@r3\nTTTTT\n+\nIIIII\n";
    // 3 reads, 8 + 3 + 5 = 16 bases.

    /// BGZF-compressed FASTQ must be detected and decoded via the parallel path,
    /// and yield exactly the same stats as the same reads read as plain text.
    #[test]
    fn bgzf_fastq_detected_and_matches_plain() {
        let dir = std::env::temp_dir();
        let bgzf_path = dir.join("wisteria_test_bgzf.fastq.gz");
        let plain_path = dir.join("wisteria_test_plain.fastq");

        // Write a real BGZF stream via noodles (no external bgzip needed).
        let mut w = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        w.write_all(FASTQ).unwrap();
        w.finish().unwrap(); // writes the BGZF EOF block

        // Plain uncompressed FASTQ baseline.
        File::create(&plain_path).unwrap().write_all(FASTQ).unwrap();

        assert!(is_bgzf(&bgzf_path).unwrap(), "bgzip output must be detected as BGZF");
        assert!(!is_bgzf(&plain_path).unwrap(), "plain FASTQ must not be BGZF");

        let (mut acc_b, stats_b) =
            process_fastq(&bgzf_path, 2, None, None, true).unwrap();
        let (mut acc_p, stats_p) =
            process_fastq(&plain_path, 2, None, None, true).unwrap();

        assert!(stats_b.bgzf_parallel, "BGZF input should decode in parallel");
        assert!(!stats_p.bgzf_parallel, "plain input should not");

        let sum_b = acc_b.calculate_summary();
        let sum_p = acc_p.calculate_summary();
        assert_eq!(sum_b.total_reads, 3);
        assert_eq!(sum_b.total_bases, 16);
        assert_eq!(sum_b.total_reads, sum_p.total_reads);
        assert_eq!(sum_b.total_bases, sum_p.total_bases);
        assert_eq!(sum_b.n50_read_length, sum_p.n50_read_length);

        let _ = std::fs::remove_file(&bgzf_path);
        let _ = std::fs::remove_file(&plain_path);
    }
}
