# Wisteria

![Wisteria Logo](resources/wisteria_logo.png)

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



## Benchmarks

Wisteria is engineered for speed, significantly outperforming standard tools on long-read datasets. Below is a benchmark comparison of running Wisteria vs. FastQC on the same dataset of **50,000 genomic long reads**:

| Tool & Input Format | Runtime (Seconds) | Speedup vs. FastQC |
| :--- | :--- | :--- |
| **FastQC** (FASTQ) | **314.49s** (5.2m) | *Baseline* |
| **Wisteria** (BAM) | **27.79s** | **11.3x faster** |
| **Wisteria** (FASTQ) | **25.37s** | **12.4x faster** |

*Note: Benchmarks were performed on 1 thread. Wisteria's parallel decompression and zero-allocation execution allow it to scale smoothly as thread count increases.*



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



## How It Works & Algorithms

Wisteria achieves its high performance through several key architectural designs and optimized algorithmic models:

### 1. Parallel BGZF Block Decompression (BAM)
BAM files are compressed using Blocked GNU Zip Format (BGZF). Decompressing these blocks is traditionally a serial CPU bottleneck. 
- **Algorithm:** Wisteria integrates `noodles::bgzf::io::MultithreadedReader`, which spins up a dedicated pool of decompression workers in the background. As the main thread requests data, BGZF blocks are fetched, sent to background threads for decompression, and returned in order. This keeps the main reader thread saturated without waiting on gzip decoding.

### 2. Rayon Batch Fold/Reduce Concurrency
Rather than processing reads sequentially or spinning up standard OS threads (which incur high context-switching and locks overhead), Wisteria utilizes **Rayon** for data-parallel task scheduling.
- **Workflow:** The main thread streams records from the decompressor, batching them in groups of 10,000. 
- **The Fold Pattern:** To avoid lock contention and memory allocation thrashing, Wisteria uses `.fold(QCAccumulator::new, ...)` to parse and accumulate statistics **thread-locally** within Rayon's worker threads. This compiles the read-level statistics (such as percentiles and density grids) locally on each CPU core.
- **The Reduce Pattern:** At the end of the batch, the thread-local accumulators are combined using a binary reduction tree (`.reduce(...)`), minimizing thread synchronization overhead and memory merges.

### 3. Single-Pass Unified Position Loop
Long reads contain tens of thousands to millions of bases per sequence. Iterating over these sequences multiple times is highly CPU-intensive due to memory cache misses and loop branching.
- **Optimization:** Wisteria replaces multiple independent sweeps (quality percentiles, absolute kb composition, GC content counts) with a **single unified loop** over the sequence. This loops through the read once, updating GC counts, 100 relative percentile bins, and absolute 1kb bins in a single pass. This dramatically increases CPU instructions-per-cycle (IPC) and speeds up the runtime by over 2x.

### 4. Zero-Allocation BAM Parsing
In standard formats, BAM sequences are stored in a 4-bit compressed form.
- **Optimization:** Rather than unpacking BAM sequences and quality scores into intermediate heap-allocated `Vec<u8>` arrays, Wisteria's custom `add_read_bam` reads and translates bases **directly from the noodles stream iterators**. This eliminates hundreds of thousands of heap allocations per file, entirely bypassing memory allocator contention.



## Output JSON Schema

The unified JSON report contains structured, downstream-ready keys:
- `summary`: High-level metrics (total/passed reads, bases, means, medians, N50, mins, maxs).
- `gc_content_distribution`: 101-element array showing read counts for each GC percentage (0% to 100%).
- `length_vs_quality_2d`: A 50x50 log-spaced length vs linear quality grid of read counts.
- `sequence_length_distribution`: An explicit 1D histogram struct containing `bin_edges` (50 bins) and the corresponding read `counts` for each length interval.
- `per_sequence_quality_distribution`: An explicit 1D histogram struct containing `bin_edges` (Q0 to Q60, 60 bins) and the corresponding read `counts` for each Phred quality score.
- `quality_by_position_percentile` / `base_content_by_position_percentile`: 100-bin arrays mapping average quality and nucleotide counts along relative read percentiles (0-100%).
- `quality_by_position_absolute` / `base_content_by_position_absolute`: 101-bin arrays mapping metrics in absolute 1kb bins (0-1kb, 1-2kb, ..., 100kb+).
