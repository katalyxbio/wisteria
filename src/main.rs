use clap::Parser;
use std::path::PathBuf;
use wisteria::debug::DebugReport;
use wisteria::io::{benchmark_decompress, process_bam, process_fastq};
use wisteria::report::generate_report;

/// Records per Rayon batch; kept in sync with the `io::process_*` readers.
const BATCH_SIZE: usize = 10_000;

/// Dispatch BAM vs FASTQ purely on the `.bam` extension.
fn is_bam_path(path: &std::path::Path) -> bool {
    path.extension()
        .map(|ext| ext.to_string_lossy().to_lowercase() == "bam")
        .unwrap_or(false)
}

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

    /// Number of worker threads for the statistics fold (Rayon pool). Defaults to all cores.
    #[arg(short, long)]
    threads: Option<usize>,

    /// Number of BGZF block-decompression workers for BAM input. Defaults to `--threads`.
    ///
    /// Reading and folding run concurrently and share the CPU, so on a busy box
    /// splitting the core budget between this and `--threads` (rather than giving
    /// both the full core count) can reduce oversubscription and speed BAM reads.
    #[arg(long)]
    decompress_threads: Option<usize>,

    /// Also render PNG charts of the report into this directory (created if
    /// absent). Rust-native rendering via `plotters`; one PNG per chart.
    #[arg(long, value_name = "DIR")]
    plots: Option<PathBuf>,

    /// Optional minimum read length filter.
    #[arg(long)]
    min_len: Option<usize>,

    /// Optional minimum average read quality score filter.
    #[arg(long)]
    min_qual: Option<f32>,

    /// Output a comprehensive debug report: per-stage timings, throughput,
    /// decompression ratio, parallelism breakdown, and memory usage.
    #[arg(long)]
    debug_stats: bool,

    /// Diagnostic (BAM only): measure raw BGZF decompression throughput without
    /// BAM record parsing, then exit. Isolates decompression cost from parsing.
    #[arg(long)]
    bench_decompress: bool,
}

fn main() {
    let args = Args::parse();
    
    let start_time = std::time::Instant::now();
    
    let threads = args.threads.unwrap_or_else(num_cpus::get);
    let decompress_threads = args.decompress_threads.unwrap_or(threads);

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
    println!("Fold threads:     {}", threads);
    if is_bam_path(&args.input) {
        println!("Decompress workers: {}", decompress_threads);
    }
    if let Some(l) = args.min_len {
        println!("Min length:       {}", l);
    }
    if let Some(q) = args.min_qual {
        println!("Min quality:      {}", q);
    }
    if let Some(dir) = &args.plots {
        println!("Plots dir:        {:?}", dir);
    }
    if args.debug_stats {
        println!("Debug Stats:      Enabled");
    }
    println!("--------------------------------");
    
    if !args.input.exists() {
        eprintln!("Error: Input file does not exist: {:?}", args.input);
        std::process::exit(1);
    }
    
    let is_bam = is_bam_path(&args.input);

    if args.bench_decompress {
        if !is_bam {
            eprintln!("Error: --bench-decompress only applies to BAM input");
            std::process::exit(1);
        }
        println!("Benchmarking raw BGZF decompression (no record parsing)...");
        match benchmark_decompress(&args.input, decompress_threads) {
            Ok((bytes, elapsed)) => {
                let secs = elapsed.as_secs_f64().max(f64::MIN_POSITIVE);
                let out_gbps = bytes as f64 / 1e9 / secs;
                let file_mb = std::fs::metadata(&args.input).map(|m| m.len()).unwrap_or(0) as f64
                    / (1024.0 * 1024.0);
                println!("Decompress workers:        {}", decompress_threads);
                println!("Uncompressed bytes:        {}", bytes);
                println!("Decompress-only wall:      {:.3?}", elapsed);
                println!("Uncompressed throughput:   {out_gbps:.2} GB/s");
                println!("Input throughput:          {:.1} MB/s", file_mb / secs);
            }
            Err(e) => {
                eprintln!("Error during decompression benchmark: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    println!("Processing starting...");

    let result = if is_bam {
        process_bam(&args.input, decompress_threads, args.min_len, args.min_qual, args.debug_stats)
    } else {
        process_fastq(&args.input, decompress_threads, args.min_len, args.min_qual, args.debug_stats)
    };

    let (mut accumulator, io_stats) = match result {
        Ok(acc) => acc,
        Err(e) => {
            eprintln!("Error processing input: {}", e);
            std::process::exit(1);
        }
    };

    let summary_start = std::time::Instant::now();
    if args.debug_stats {
        println!("Generating statistics summary...");
    }
    let summary = accumulator.calculate_summary();
    let summary_duration = summary_start.elapsed();

    // Retained per-read vectors drive the accumulator's live memory; capture
    // their length before `generate_report` consumes the accumulator.
    let retained_reads = accumulator.read_lengths.len();

    let report_start = std::time::Instant::now();
    println!("Writing unified JSON report...");
    if args.plots.is_some() {
        println!("Rendering plots...");
    }
    let plots = match generate_report(accumulator, &args.output, args.plots.as_deref()) {
        Ok(plots) => plots,
        Err(e) => {
            eprintln!("Error writing report: {}", e);
            std::process::exit(1);
        }
    };
    let report_duration = report_start.elapsed();

    if let Some(dir) = &args.plots {
        println!("Wrote {} plot(s) to {:?}:", plots.len(), dir);
        for name in &plots {
            println!("  - {}", name);
        }
    }

    let duration = start_time.elapsed();

    if args.debug_stats {
        let report = DebugReport {
            input: &args.input,
            output: &args.output,
            is_bam,
            threads,
            decompress_threads,
            min_len: args.min_len,
            min_qual: args.min_qual,
            batch_size: BATCH_SIZE,
            io: io_stats,
            summary_time: summary_duration,
            report_time: report_duration,
            total_time: duration,
            input_file_bytes: std::fs::metadata(&args.input).ok().map(|m| m.len()),
            total_reads: summary.total_reads,
            total_bases: summary.total_bases,
            passed_reads: summary.passed_filters_reads,
            passed_bases: summary.passed_filters_bases,
            retained_reads,
        };
        print!("{}", report.render());
    }

    println!("\nQC Summary Metrics:");
    println!("==================================================");
    println!("Total Reads:               {}", summary.total_reads);
    println!("Total Bases:               {}", summary.total_bases);
    if summary.skipped_alignments > 0 {
        println!("Skipped (secondary/suppl.):{}", summary.skipped_alignments);
    }
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
