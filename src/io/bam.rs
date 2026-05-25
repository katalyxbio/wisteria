use std::fs::File;
use std::path::Path;
use std::num::NonZeroUsize;
use std::io::Write;
use std::time::Duration;
use crate::stats::QCAccumulator;
use noodles::bam;
use noodles::bgzf;
use rayon::prelude::*;

pub fn process_bam(
    path: &Path,
    num_threads: usize,
    min_len: Option<usize>,
    min_qual: Option<f32>,
    debug: bool,
) -> Result<QCAccumulator, Box<dyn std::error::Error>> {
    let start_time = std::time::Instant::now();
    let file = File::open(path)?;
    
    // Setup multi-threaded BGZF decompression
    let worker_count = NonZeroUsize::new(num_threads).unwrap_or(NonZeroUsize::new(1).unwrap());
    let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    
    let _header = reader.read_header()?;
    
    let mut accumulator = QCAccumulator::new();
    let batch_size = 10_000;
    let mut batch = Vec::with_capacity(batch_size);
    let mut total_processed = 0;
    
    let mut total_io_time = Duration::default();
    let mut total_proc_time = Duration::default();
    let mut total_merge_time = Duration::default();
    
    let mut iter = reader.records();
    
    loop {
        let io_start = std::time::Instant::now();
        let mut batch_filled = false;
        
        for _ in 0..batch_size {
            if let Some(result) = iter.next() {
                let record = result?;
                batch.push(record);
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
        let local_acc = process_batch(&batch, min_len, min_qual);
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
        println!("\n[DEBUG BAM TIMINGS]");
        println!("--------------------------------------------------");
        println!("BGZF I/O & Parse:          {:.3?}", total_io_time);
        println!("Rayon Concurrency (Fold):  {:.3?}", total_proc_time);
        println!("Main Thread Accumulation:  {:.3?}", total_merge_time);
        println!("Total File Time:           {:.3?}", start_time.elapsed());
        println!("--------------------------------------------------");
    }
    
    Ok(accumulator)
}

fn process_batch(
    records: &[bam::Record],
    min_len: Option<usize>,
    min_qual: Option<f32>,
) -> QCAccumulator {
    records
        .par_iter()
        .fold(QCAccumulator::new, |mut acc, record| {
            acc.add_read_bam(record, min_len, min_qual);
            acc
        })
        .reduce(QCAccumulator::new, |mut a, b| {
            a.merge(b);
            a
        })
}
