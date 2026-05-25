use clap::Parser;
use std::path::PathBuf;
use wisteria::io::{process_bam, process_fastq};
use wisteria::report::generate_report;

#[derive(Parser, Debug)]
#[command(name = "Wisteria")]
#[command(author = "DeepMind Pair Programmer")]
#[command(version = "0.1.0")]
#[command(about = "Ultra-fast genomic long-read quality check", long_about = None)]
struct Args {
    /// Path to input FASTQ or BAM file. FASTQ files can be gzipped (.gz).
    #[arg(short, long)]
    input: PathBuf,

    /// Path to output unified JSON report.
    #[arg(short, long, default_value = "wisteria_report.json")]
    output: PathBuf,

    /// Number of worker threads for parallel processing (BAM decompression and statistics compilation).
    #[arg(short, long)]
    threads: Option<usize>,

    /// Optional minimum read length filter.
    #[arg(long)]
    min_len: Option<usize>,

    /// Optional minimum average read quality score filter.
    #[arg(long)]
    min_qual: Option<f32>,

    /// Track step timings and output comprehensive debug statistics.
    #[arg(long)]
    debug: bool,
}

fn main() {
    let args = Args::parse();
    
    let start_time = std::time::Instant::now();
    
    let threads = args.threads.unwrap_or_else(num_cpus::get);
    
    // Set up global Rayon thread pool if requested threads != default
    if args.threads.is_some() {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
    
    println!("Wisteria long-read Quality Check");
    println!("--------------------------------");
    println!("Input file:       {:?}", args.input);
    println!("Output report:    {:?}", args.output);
    println!("Threads:          {}", threads);
    if let Some(l) = args.min_len {
        println!("Min length:       {}", l);
    }
    if let Some(q) = args.min_qual {
        println!("Min quality:      {}", q);
    }
    if args.debug {
        println!("Debug Mode:       Enabled");
    }
    println!("--------------------------------");
    
    if !args.input.exists() {
        eprintln!("Error: Input file does not exist: {:?}", args.input);
        std::process::exit(1);
    }
    
    let is_bam = args.input.extension()
        .map(|ext| ext.to_string_lossy().to_lowercase() == "bam")
        .unwrap_or(false);
        
    println!("Processing starting...");
    
    let proc_start = std::time::Instant::now();
    let result = if is_bam {
        process_bam(&args.input, threads, args.min_len, args.min_qual, args.debug)
    } else {
        process_fastq(&args.input, args.min_len, args.min_qual, args.debug)
    };
    
    let mut accumulator = match result {
        Ok(acc) => acc,
        Err(e) => {
            eprintln!("Error processing input: {}", e);
            std::process::exit(1);
        }
    };
    let proc_duration = proc_start.elapsed();
    
    let summary_start = std::time::Instant::now();
    if args.debug {
        println!("Generating statistics summary...");
    }
    let summary = accumulator.calculate_summary();
    let summary_duration = summary_start.elapsed();
    
    let report_start = std::time::Instant::now();
    if args.debug {
        println!("Writing unified JSON report...");
    } else {
        println!("Writing unified JSON report...");
    }
    if let Err(e) = generate_report(accumulator, &args.output) {
        eprintln!("Error writing report: {}", e);
        std::process::exit(1);
    }
    let report_duration = report_start.elapsed();
    
    let duration = start_time.elapsed();
    
    if args.debug {
        println!("\n[DEBUG EXECUTIVE TIMINGS]");
        println!("--------------------------------------------------");
        println!("File Processing:           {:.3?}", proc_duration);
        println!("Summary Stats Math:        {:.3?}", summary_duration);
        println!("JSON Serialization & I/O:  {:.3?}", report_duration);
        println!("Grand Total Execution:     {:.3?}", duration);
        println!("--------------------------------------------------");
    }
    
    println!("\nQC Summary Metrics:");
    println!("==================================================");
    println!("Total Reads:               {}", summary.total_reads);
    println!("Total Bases:               {}", summary.total_bases);
    println!("Passed Filters Reads:      {}", summary.passed_filters_reads);
    if summary.total_reads > 0 {
        let pct = (summary.passed_filters_reads as f64 / summary.total_reads as f64) * 100.0;
        println!("Passed Filters Reads (%):  {:.2}%", pct);
    }
    println!("Passed Filters Bases:      {}", summary.passed_filters_bases);
    println!("Mean Read Length:          {:.1} bp", summary.mean_read_length);
    println!("Median Read Length:        {} bp", summary.median_read_length);
    println!("N50 Read Length:           {} bp", summary.n50_read_length);
    println!("Min Read Length:           {} bp", summary.min_read_length);
    println!("Max Read Length:           {} bp", summary.max_read_length);
    println!("Mean Read Quality:         {:.2}", summary.mean_read_quality);
    println!("Median Read Quality:       {:.2}", summary.median_read_quality);
    println!("Min Read Quality:          {:.2}", summary.min_read_quality);
    println!("Max Read Quality:          {:.2}", summary.max_read_quality);
    println!("==================================================");
    println!("Completed in {:.2?}", duration);
}
