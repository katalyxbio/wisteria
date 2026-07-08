use std::fs::File;
use std::path::Path;
use std::num::NonZeroUsize;
use std::io::{BufRead, Write};
use std::sync::mpsc::sync_channel;
use std::time::{Duration, Instant};
use crate::debug::IoStats;
use crate::stats::QCAccumulator;
use noodles::bam;
use noodles::bgzf;
use rayon::prelude::*;

const BATCH_SIZE: usize = 10_000;
/// Also cap a batch by bytes so a run of very long reads can't blow up memory.
const MAX_BATCH_BYTES: usize = 128 << 20; // 128 MiB
/// Batches the reader may run ahead of the folder; bounds memory and provides
/// backpressure so we never buffer the whole file.
const CHANNEL_DEPTH: usize = 2;

/// A batch of raw BAM record bytes plus the `(offset, len)` span of each record
/// within the buffer. Records are parsed zero-copy via `bam::RecordRef` in the
/// fold, keeping the single-threaded producer's work to bulk reads + framing.
type BamBatch = (Vec<u8>, Vec<(usize, usize)>);

pub fn process_bam(
    path: &Path,
    num_threads: usize,
    min_len: Option<usize>,
    min_qual: Option<f32>,
    debug_stats: bool,
) -> Result<(QCAccumulator, IoStats), Box<dyn std::error::Error>> {
    let file = File::open(path)?;

    // Setup multi-threaded BGZF decompression
    let worker_count = NonZeroUsize::new(num_threads).unwrap_or(NonZeroUsize::new(1).unwrap());
    let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);

    let _header = reader.read_header()?;

    // Drop down to the raw BGZF byte stream (positioned at the first record) so
    // the producer can frame records itself instead of paying the per-record
    // allocate + zero-fill + copy cost of the `records()` iterator.
    let mut inner = reader.into_inner();

    // Pipeline: this thread frames raw record spans while a consumer thread
    // parses (zero-copy) and folds them, overlapping I/O with compute.
    let (tx, rx) = sync_channel::<BamBatch>(CHANNEL_DEPTH);

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
            let (buf, spans) = match rx.recv() {
                Ok(batch) => batch,
                Err(_) => break,
            };
            recv_wait += recv_start.elapsed();

            let proc_start = Instant::now();
            let local_acc = process_batch(&buf, &spans, min_len, min_qual);
            proc_time += proc_start.elapsed();

            let merge_start = Instant::now();
            total_processed += spans.len();
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

    // Producer: frame raw record byte-spans on this thread and hand them off.
    let mut total_io_time = Duration::default();
    let mut send_wait = Duration::default();
    let mut read_error: Option<Box<dyn std::error::Error>> = None;
    let mut carry: Vec<u8> = Vec::new(); // leftover bytes straddling a batch boundary

    'outer: loop {
        let io_start = Instant::now();

        let mut buf = std::mem::take(&mut carry);
        let mut spans: Vec<(usize, usize)> = Vec::with_capacity(BATCH_SIZE);
        let mut pos = 0usize;
        let mut eof = false;

        loop {
            // Frame as many complete records as are already buffered.
            while spans.len() < BATCH_SIZE && pos < MAX_BATCH_BYTES {
                if buf.len() - pos < 4 {
                    break;
                }
                let block_size = u32::from_le_bytes([
                    buf[pos],
                    buf[pos + 1],
                    buf[pos + 2],
                    buf[pos + 3],
                ]) as usize;
                if buf.len() - pos - 4 < block_size {
                    break; // record not fully buffered yet
                }
                spans.push((pos + 4, block_size));
                pos += 4 + block_size;
            }

            if spans.len() >= BATCH_SIZE || pos >= MAX_BATCH_BYTES {
                break;
            }

            // Need more bytes: append the next decompressed block (single copy).
            let src = match inner.fill_buf() {
                Ok(src) => src,
                Err(e) => {
                    read_error = Some(Box::new(e));
                    eof = true;
                    break;
                }
            };
            if src.is_empty() {
                eof = true;
                break;
            }
            let n = src.len();
            buf.extend_from_slice(src);
            inner.consume(n);
        }

        // Bytes after the last framed record carry over to the next batch.
        if pos < buf.len() {
            carry = buf[pos..].to_vec();
            buf.truncate(pos);
        }

        total_io_time += io_start.elapsed();

        if !spans.is_empty() {
            let send_start = Instant::now();
            let sent = tx.send((buf, spans));
            send_wait += send_start.elapsed();
            if sent.is_err() {
                break 'outer; // consumer hung up
            }
        }
        if eof {
            break;
        }
    }
    drop(tx); // close the channel so the consumer's loop terminates

    let (accumulator, proc_time, merge_time, recv_wait, batches) = consumer
        .join()
        .map_err(|_| Box::<dyn std::error::Error>::from("BAM consumer thread panicked"))?;

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
        bgzf_parallel: true, // BAM is always BGZF, decoded in parallel
    };

    Ok((accumulator, stats))
}

fn process_batch(
    buf: &[u8],
    spans: &[(usize, usize)],
    min_len: Option<usize>,
    min_qual: Option<f32>,
) -> QCAccumulator {
    spans
        .par_iter()
        .fold(QCAccumulator::new, |mut acc, &(offset, len)| {
            // Zero-copy view over the record bytes; skip anything too short to
            // be a valid record (RecordRef::new requires >= 32 bytes).
            if let Some(record) = bam::RecordRef::new(&buf[offset..offset + len]) {
                acc.add_read_bam(&record, min_len, min_qual);
            }
            acc
        })
        .reduce(QCAccumulator::new, |mut a, b| {
            a.merge(b);
            a
        })
}

/// Diagnostic: drain the raw BGZF stream (decompress only, no BAM record
/// parsing) to isolate pure decompression throughput from record parsing.
/// Returns the number of uncompressed bytes read and the elapsed wall time.
pub fn benchmark_decompress(
    path: &Path,
    num_threads: usize,
) -> Result<(u64, Duration), Box<dyn std::error::Error>> {
    use std::io::Read;

    let file = File::open(path)?;
    let worker_count = NonZeroUsize::new(num_threads).unwrap_or(NonZeroUsize::new(1).unwrap());
    let mut reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);

    let start = Instant::now();
    let mut buf = vec![0u8; 1 << 20];
    let mut total: u64 = 0;
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        total += n as u64;
    }
    Ok((total, start.elapsed()))
}
