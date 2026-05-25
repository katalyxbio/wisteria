# Wisteria

Wisteria is a Rust-based, ultra-fast quality control (QC) tool designed specifically for genomic long-read datasets (such as Oxford Nanopore and PacBio). It takes both BAM and FASTQ (raw or gzipped) as input, handles multi-threading elegantly, and outputs a single unified JSON report containing comprehensive statistics and binned distributions.

## Key Features

- **Blazingly Fast & Scalable:** Leverages a data-parallel MapReduce architecture powered by **Rayon** for thread-safe accumulator reductions.
- **Advanced BAM Concurrency:** Utilizes the pure-Rust **`noodles`** library with `bgzf::io::MultithreadedReader` to handle multi-threaded block decompression.
- **Robust FASTQ Parsing:** Transparently detects and decompresses Gzip files (.gz) using zero-copy stream processing via **`needletail`**.
- **Long-Read Tailored Statistics:**
  - Standard metrics adapted to long reads: N50 read length, exact median length/quality, mean, min, and max values.
  - **Percentile Binning:** Dividers that split each read into 100 sections to track quality decay and base composition along the relative length.
  - **Absolute Binning:** Tracks quality and base composition in 1kb intervals up to 100kb+ for direct sequence decay observation.
  - **2D Density Matrix:** Generates a log-spaced length vs. quality density matrix to capture exact sequence profile clusters (ideal for 2D heatmaps).
- **Optional Read Filtering:** Adds optional, on-the-fly pre-filtering with `--min-len` and `--min-qual` to clean up datasets before calculating statistics.
- **Dynamic Progress Counter:** Displays a real-time processed reads progress counter in the console during execution.

---

## Installation

### Prerequisites
Make sure you have Rust and Cargo installed. If not, install them via [rustup.rs](https://rustup.rs/):
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Global Installation
Install `wisteria` directly to your local Cargo bin directory (`~/.cargo/bin`):
```bash
# From the wisteria workspace directory:
cargo install --path .
```
Verify the installation was successful and is available globally in your PATH:
```bash
wisteria --help
```

---

## Usage

```text
Ultra-fast genomic long-read quality check

Usage: wisteria [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>        Path to input FASTQ or BAM file. FASTQ files can be gzipped (.gz)
  -o, --output <OUTPUT>      Path to output unified JSON report [default: wisteria_report.json]
  -t, --threads <THREADS>    Number of worker threads for parallel processing (BAM decompression and statistics compilation)
      --min-len <MIN_LEN>    Optional minimum read length filter
      --min-qual <MIN_QUAL>  Optional minimum average read quality score filter
  -h, --help                 Print help
  -V, --version              Print version
```

### Running Examples

1. **Analyze a gzipped FASTQ file:**
   ```bash
   wisteria --input sample.fastq.gz --output sample_qc.json
   ```

2. **Analyze a BAM file utilizing 8 decompression and worker threads:**
   ```bash
   wisteria --input aligned_reads.bam --output bam_qc.json --threads 8
   ```

3. **Apply quality and length filters during processing:**
   ```bash
   # Skips reads shorter than 1,000 bp or with an average Phred quality score below Q9
   wisteria --input nanopore.fastq.gz --min-len 1000 --min-qual 9.0
   ```

---

## Output JSON Schema

The unified JSON report contains structured, downstream-ready keys:
- `summary`: High-level metrics (total/passed reads, bases, means, medians, N50, mins, maxs).
- `gc_content_distribution`: 101-element array showing read counts for each GC percentage (0% to 100%).
- `length_vs_quality_2d`: A 50x50 log-spaced length vs linear quality grid of read counts.
- `quality_by_position_percentile` / `base_content_by_position_percentile`: 100-bin arrays mapping average quality and nucleotide counts along relative read percentiles (0-100%).
- `quality_by_position_absolute` / `base_content_by_position_absolute`: 101-bin arrays mapping metrics in absolute 1kb bins (0-1kb, 1-2kb, ..., 100kb+).
