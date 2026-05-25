use std::path::Path;
use std::time::Duration;
use crate::stats::QCAccumulator;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::io::Write;

pub fn process_fastq(
    path: &Path,
    min_len: Option<usize>,
    min_qual: Option<f32>,
    debug: bool,
) -> Result<QCAccumulator, Box<dyn std::error::Error>> {
    let start_time = std::time::Instant::now();
    let mut reader = parse_fastx_file(path)?;
    
    let mut accumulator = QCAccumulator::new();
    let batch_size = 10_000;
    let mut batch = Vec::with_capacity(batch_size);
    let mut total_processed = 0;
    
    let mut total_io_time = Duration::default();
    let mut total_proc_time = Duration::default();
    let mut total_merge_time = Duration::default();
    
    loop {
        let io_start = std::time::Instant::now();
        let mut batch_filled = false;
        
        for _ in 0..batch_size {
            if let Some(record) = reader.next() {
                let seqrec = record?;
                let seq = seqrec.seq().to_vec();
                let qual = seqrec.qual().map(|q| q.to_vec()).unwrap_or_else(Vec::new);
                batch.push((seq, qual));
            } else {
                batch_filled = true;
                break;
            }
        }
        total_io_time += io_start.elapsed();
        
        if batch.is_empty() {
            break;
        }
        
        let proc_start = std::time::Instant::now();
        let local_acc = process_fastq_batch(&batch, min_len, min_qual);
        total_proc_time += proc_start.elapsed();
        
        let merge_start = std::time::Instant::now();
        let count = batch.len();
        accumulator.merge(local_acc);
        total_processed += count;
        total_merge_time += merge_start.elapsed();
        
        if !debug {
            print!("\rProcessed {} reads...", total_processed);
            let _ = std::io::stdout().flush();
        }
        
        batch.clear();
        
        if batch_filled {
            break;
        }
    }
    
    if !debug && total_processed > 0 {
        println!("\rProcessed {} reads. Done!", total_processed);
    }
    
    if debug {
        println!("\n[DEBUG FASTQ TIMINGS]");
        println!("--------------------------------------------------");
        println!("Decompress & Parse:        {:.3?}", total_io_time);
        println!("Rayon Concurrency (Fold):  {:.3?}", total_proc_time);
        println!("Main Thread Accumulation:  {:.3?}", total_merge_time);
        println!("Total File Time:           {:.3?}", start_time.elapsed());
        println!("--------------------------------------------------");
    }
    
    Ok(accumulator)
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
