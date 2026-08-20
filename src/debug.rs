//! Debug/profiling instrumentation for `--debug-stats`.
//!
//! Collects per-stage timings, throughput, and memory figures so the operator
//! can see exactly where wall-clock and allocations go and which stages are
//! worth optimizing (e.g. the serial merge, or the retained per-read vectors).

use std::fmt::Write as _;
use std::path::Path;
use std::time::Duration;

/// Timings and counters gathered while reading and processing an input file.
///
/// Populated by the `io::process_*` functions and consumed by
/// [`DebugReport`] to render the final `--debug-stats` block.
#[derive(Debug, Default, Clone)]
pub struct IoStats {
    /// Time spent pulling and parsing records from the (decompressed) stream.
    pub io_time: Duration,
    /// Time spent in the Rayon fold/reduce over each batch.
    pub proc_time: Duration,
    /// Time spent merging each batch's accumulator into the running total (serial).
    pub merge_time: Duration,
    /// Number of batches processed.
    pub batches: u64,
    /// Producer time blocked handing batches to the consumer (channel full →
    /// the folder is the bottleneck).
    pub send_wait: Duration,
    /// Consumer time blocked waiting for the next batch (channel empty → the
    /// reader/producer is the bottleneck).
    pub recv_wait: Duration,
    /// Whether input decompression ran in parallel (BGZF block decode). True for
    /// BAM and for BGZF-compressed FASTQ; false for plain single-stream gzip or
    /// uncompressed input, which decode serially.
    pub bgzf_parallel: bool,
}

/// Peak resident set size in kB (`VmHWM`), read from `/proc/self/status`.
///
/// Returns `None` on non-Linux platforms or if the field is unavailable.
pub fn peak_rss_kb() -> Option<u64> {
    read_status_field("VmHWM:")
}

/// Current resident set size in kB (`VmRSS`), read from `/proc/self/status`.
///
/// Returns `None` on non-Linux platforms or if the field is unavailable.
pub fn current_rss_kb() -> Option<u64> {
    read_status_field("VmRSS:")
}

fn read_status_field(field: &str) -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix(field) {
            // Format: "VmHWM:\t  412345 kB"
            return rest.split_whitespace().next()?.parse().ok();
        }
    }
    None
}

/// Everything needed to render the comprehensive `--debug-stats` report.
pub struct DebugReport<'a> {
    pub input: &'a Path,
    pub output: &'a Path,
    pub is_bam: bool,
    pub threads: usize,
    pub decompress_threads: usize,
    pub min_len: Option<usize>,
    pub min_qual: Option<f32>,
    pub batch_size: usize,

    pub io: IoStats,
    pub summary_time: Duration,
    pub report_time: Duration,
    pub total_time: Duration,

    pub input_file_bytes: Option<u64>,
    pub total_reads: u64,
    pub total_bases: u64,
    pub passed_reads: u64,
    pub passed_bases: u64,
    /// Number of reads retained in the per-read `read_lengths`/`average_qualities`
    /// vectors (== passed reads); drives the accumulator's live memory.
    pub retained_reads: usize,
}

impl DebugReport<'_> {
    /// Format the full report as a printable string.
    pub fn render(&self) -> String {
        let mut o = String::new();
        let total_secs = self.total_time.as_secs_f64().max(f64::MIN_POSITIVE);

        let _ = writeln!(o, "\n[DEBUG STATS]");
        let _ = writeln!(o, "==================================================");

        let _ = writeln!(o, "-- Run Configuration --");
        let _ = writeln!(o, "Format:               {}", if self.is_bam { "BAM" } else { "FASTQ" });
        let _ = writeln!(o, "Input:                {:?}", self.input);
        let _ = writeln!(o, "Output:               {:?}", self.output);
        let _ = writeln!(o, "Fold threads:         {}", self.threads);
        if self.is_bam || self.io.bgzf_parallel {
            let _ = writeln!(o, "Decompress workers:   {}", self.decompress_threads);
        }
        let _ = writeln!(o, "Batch size:           {}", self.batch_size);
        let _ = writeln!(o, "Min length filter:    {}", opt(self.min_len));
        let _ = writeln!(o, "Min quality filter:   {}", opt(self.min_qual));

        let _ = writeln!(o, "\n-- Input & Decompression --");
        let decode_mode = if self.io.bgzf_parallel {
            "parallel BGZF block decode"
        } else {
            "serial gzip / uncompressed"
        };
        let _ = writeln!(o, "Decompress mode:      {decode_mode}");
        match self.input_file_bytes {
            Some(bytes) => {
                let _ = writeln!(o, "Input file size:      {}", human_bytes(bytes));
                if bytes > 0 {
                    let ratio = self.total_bases as f64 / bytes as f64;
                    let _ = writeln!(o, "Uncompressed bases:   {}", self.total_bases);
                    let _ = writeln!(o, "Bases / input byte:   {ratio:.2}x");
                    let mb = bytes as f64 / (1024.0 * 1024.0);
                    let _ = writeln!(o, "Input throughput:     {:.1} MB/s", mb / total_secs);
                }
            }
            None => {
                let _ = writeln!(o, "Input file size:      n/a");
            }
        }

        let _ = writeln!(o, "\n-- Throughput & Filtering --");
        let _ = writeln!(o, "Reads/sec:            {:.0}", self.total_reads as f64 / total_secs);
        let _ = writeln!(o, "Bases/sec:            {:.0}", self.total_bases as f64 / total_secs);
        let _ = writeln!(o, "Total reads:          {}", self.total_reads);
        let _ = writeln!(o, "Passed filters:       {} reads / {} bases", self.passed_reads, self.passed_bases);
        if self.total_reads > 0 {
            let pass_pct = self.passed_reads as f64 / self.total_reads as f64 * 100.0;
            let _ = writeln!(o, "Pass rate:            {pass_pct:.2}%");
        }

        let _ = writeln!(o, "\n-- Stage Timings --");
        let _ = writeln!(o, "{:<28} {:>10}  {:>7}", "stage", "time", "% total");
        self.stage_line(&mut o, "I/O & parse", self.io.io_time, total_secs);
        self.stage_line(&mut o, "Rayon fold/reduce", self.io.proc_time, total_secs);
        self.stage_line(&mut o, "Batch merge (serial)", self.io.merge_time, total_secs);
        self.stage_line(&mut o, "Summary math", self.summary_time, total_secs);
        self.stage_line(&mut o, "JSON serialize & write", self.report_time, total_secs);
        let _ = writeln!(o, "{:-<48}", "");
        self.stage_line(&mut o, "Grand total", self.total_time, total_secs);
        let _ = writeln!(o, "(I/O & fold run concurrently; their times overlap and sum to more than the grand total.)");

        let _ = writeln!(o, "\n-- Parallelism --");
        let _ = writeln!(o, "Batches processed:    {}", self.io.batches);
        if self.io.batches > 0 {
            let avg = self.io.proc_time.as_secs_f64() / self.io.batches as f64;
            let _ = writeln!(o, "Avg fold time/batch:  {:.3} ms", avg * 1000.0);
        }
        let par = self.io.proc_time + self.io.merge_time;
        if !par.is_zero() {
            // The merge is serial main-thread work; a high share caps scaling.
            let merge_share = self.io.merge_time.as_secs_f64() / par.as_secs_f64() * 100.0;
            let _ = writeln!(o, "Serial merge share:   {merge_share:.1}% of fold+merge");
        }
        // Pipeline balance: which side stalls waiting for the other.
        let _ = writeln!(
            o,
            "Consumer recv-stall:  {:.3?}  (folder idle, waiting on the reader)",
            self.io.recv_wait
        );
        let _ = writeln!(
            o,
            "Producer send-stall:  {:.3?}  (reader idle, waiting on the folder)",
            self.io.send_wait
        );

        let _ = writeln!(o, "\n-- Memory --");
        let _ = writeln!(o, "Peak RSS (VmHWM):     {}", opt_kb(peak_rss_kb()));
        let _ = writeln!(o, "Current RSS (VmRSS):  {}", opt_kb(current_rss_kb()));
        // read_lengths: Vec<u32> + average_qualities: Vec<f32> = 8 bytes/read retained.
        let retained_bytes = self.retained_reads as u64 * 8;
        let _ = writeln!(
            o,
            "Retained per-read:    {} ({} reads x 8 B: read_lengths + average_qualities)",
            human_bytes(retained_bytes),
            self.retained_reads
        );

        let _ = writeln!(o, "==================================================");
        o
    }

    fn stage_line(&self, o: &mut String, name: &str, d: Duration, total_secs: f64) {
        let pct = d.as_secs_f64() / total_secs * 100.0;
        let _ = writeln!(o, "{name:<28} {:>10.3?}  {pct:>6.1}%", d);
    }
}

fn opt<T: std::fmt::Display>(v: Option<T>) -> String {
    match v {
        Some(v) => v.to_string(),
        None => "none".to_string(),
    }
}

fn opt_kb(kb: Option<u64>) -> String {
    match kb {
        Some(kb) => human_bytes(kb * 1024),
        None => "n/a".to_string(),
    }
}

fn human_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KB", "MB", "GB", "TB"];
    let mut val = bytes as f64;
    let mut unit = 0;
    while val >= 1024.0 && unit < UNITS.len() - 1 {
        val /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{val:.2} {}", UNITS[unit])
    }
}
