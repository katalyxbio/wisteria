# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Wisteria is an ultra-fast quality-control CLI for genomic **long-read** datasets (Oxford Nanopore, PacBio). It accepts BAM and FASTQ (raw or gzipped) input and emits a single unified JSON report of summary stats and binned distributions. Rust 2024 edition.

## Commands

```bash
cargo build --release          # optimized binary (target/release/wisteria)
cargo run -- --input <file>    # run in debug build
cargo test                     # run unit tests (see src/stats/accumulator.rs)
cargo test test_n50_and_median # run a single test by name
cargo install --path .         # install to ~/.cargo/bin
```

Run the tool: `wisteria --input reads.fastq.gz --output report.json --threads 8`.
Add `--debug-stats` for a comprehensive profiling report: per-stage timings (I/O vs. Rayon fold vs. serial merge), throughput, decompression ratio, parallelism breakdown, and memory (peak/current RSS from `/proc`). See `src/debug.rs`.

The `release` profile uses fat LTO, `codegen-units = 1`, and `panic = "abort"` — release builds are slow; prefer `cargo build` / `cargo test` while iterating.

## Architecture

Pipeline: **read → batch → parallel fold/reduce into an accumulator → summarize → serialize JSON.**

- `src/main.rs` — CLI parsing (clap), dispatches to BAM vs. FASTQ by file extension (`.bam` → BAM, everything else → FASTQ), prints summary, writes report.
- `src/io/` — input readers, one per format, each returning a fully-merged `QCAccumulator`:
  - `bam.rs`: reads the header via the `noodles` BAM reader over `bgzf::io::MultithreadedReader` (parallel block decompression), then drops to the raw BGZF byte stream via `into_inner()`. The producer **frames raw record byte-spans itself** — copying each decompressed block once into a batch buffer via `BufRead::fill_buf` and recording each record's `(offset, len)` — instead of the `records()` iterator, which allocates + zero-fills + copies a `Vec` per record (the old single-threaded bottleneck). Records are parsed **zero-copy** in the fold via `bam::RecordRef::new(&buf[span])`. Batches are capped at 10,000 records or `MAX_BATCH_BYTES` (128 MiB), whichever comes first.
  - `fastq.rs`: sniffs the first 18 bytes for the BGZF signature (`is_bgzf`: gzip magic + `FEXTRA` + `BC` subfield). **BGZF-compressed FASTQ** (from `bgzip`) is decoded in parallel via the same `bgzf::io::MultithreadedReader` as the BAM path, then handed to `needletail::parse_fastx_reader` (which parses the already-decompressed plaintext and adds no gzip layer). **Plain single-stream `.fastq.gz`** falls back to `needletail::parse_fastx_file` — a *serial* DEFLATE decode that is the sole bottleneck on large gzip inputs (measured: ~525s decode vs ~35s fold on a 50 GB file; parallel BGZF removes the wall). Batches of 10,000 owned `(seq, qual)` pairs either way. `--decompress-threads` sizes the BGZF worker pool (ignored for plain gzip, which can't parallelize).
  - **Both are pipelined** (producer/consumer): the calling thread reads/frames batches and sends them over a bounded `sync_channel` (depth 2) to a spawned consumer thread that folds and merges — so decompression/framing overlaps with the statistics fold. The consumer owns the running `QCAccumulator` and returns it (plus timings) on join. (The FASTQ reader is `!Send` so it must stay on the producer thread; the BAM producer keeps the BGZF reader there too.)
  - Each batch is folded with Rayon `par_iter().fold(QCAccumulator::new, ...).reduce(...)` — statistics accumulate **thread-locally** per core, then combine via a binary reduction tree. The consumer merges each batch's result serially (cheap: measured <0.1% of runtime). Profiling showed the BAM wall is producer-bound on record reading, not decompression (`--bench-decompress` isolates the two), which is why the producer avoids per-record work and defers parsing into the parallel fold.
- `src/stats/accumulator.rs` — `QCAccumulator` is the core data structure and the unit of the MapReduce. `add_read` (FASTQ) and `add_read_bam` (BAM) each do a **single unified pass** over the sequence, updating GC counts, 100 relative-percentile bins, and absolute 1 kb bins together. `merge()` combines two accumulators (required by Rayon `reduce`). `calculate_summary()` computes mean/median/N50/min/max (mutates: sorts `read_lengths` and `average_qualities` in place).
- `src/stats/binned_stats.rs` — `PositionStats` (per-bin base + quality counters), `BinnedStatsPercentile` (100 relative bins), `BinnedStatsAbsolute` (1 kb absolute bins up to 100 kb+).
- `src/report/json.rs` — builds `UnifiedQCReport`, averages the binned sums into per-position vectors, computes the 1D length/quality histograms, and writes pretty JSON. See README "Output JSON Schema" for the full key list. `generate_report(acc, output_path, plots_dir: Option<&Path>)` returns the list of plot filenames written (empty when `plots_dir` is `None`); when `Some`, it renders after the JSON via `plots::generate_plots`.
- `src/report/plots.rs` — `--plots <dir>` PNG rendering via **`plotters`** (pure-Rust: `bitmap_backend` + `bitmap_encoder` + `ab_glyph`; no freetype/font-kit C deps). Renders 7 charts from `&UnifiedQCReport` (length hist, quality hist, GC distribution, quality-by-position percentile + absolute, base-content-by-position, and the length-vs-quality heatmap). `ab_glyph` ships **no default font**, so a Roboto TTF (`resources/Roboto-Regular.ttf`, Apache-2.0) is `include_bytes!`-embedded and registered once via `register_font` under the `sans-serif`/`serif`/`monospace` families a `std::sync::Once`; without this, rendering fails with `FontUnavailable`.
- `src/debug.rs` — `--debug-stats` instrumentation: `IoStats` (per-stage timings, batch count, and producer `send_wait` / consumer `recv_wait` pipeline-stall timings, populated by the `io::process_*` readers and returned alongside the accumulator), `/proc/self/status` RSS helpers (Linux-only, `n/a` elsewhere), and `DebugReport::render` for the formatted report. `recv_wait >> send_wait` means the reader is the bottleneck; the reverse means the fold is. `io::benchmark_decompress` (behind `--bench-decompress`) drains the raw BGZF stream with no record parsing to isolate decompression cost from BAM parsing.

## Key conventions & gotchas

- **Quality-score encoding differs by format.** FASTQ quality bytes are ASCII Phred+33 (subtract 33), BAM stores raw Phred values. `add_read` takes an `is_fastq_qual: bool` flag; `add_read_bam` always treats scores as raw. Preserve this when touching either path.
- **`--min-len`/`--min-qual` filtering is inline and asymmetric.** In both `add_read` and `add_read_bam`, `total_reads`/`total_bases` are incremented *before* the filter check, but the early `return` on failure happens *before* any read is pushed to `read_lengths`/`average_qualities` or folded into the GC/positional/2D bins. So `total_reads` counts every read seen (minus BAM secondary/supplementary), while every other stat (medians, N50, histograms, per-position bins) reflects only reads that passed both filters — mirrored by `passed_filters_reads`/`passed_filters_bases`.
- **`add_read_bam` skips secondary (`0x100`) and supplementary (`0x800`) alignments** (checked via `record.flags()` before any counting) — QC describes reads, not alignments, and those records carry no/partial `SEQ`. They're tallied in `skipped_alignments` (reported in the summary) but excluded from `total_reads`/stats. Unmapped primary reads (`0x4`) are kept — the common case for QC on unaligned long reads. FASTQ has no such concept.
- **`add_read` (takes `&[u8]` seq/qual from FASTQ) and `add_read_bam` (takes a `&bam::RecordRef`) are near-duplicate loops** that must stay in sync — the BAM version reads directly from the record's zero-copy sequence/quality iterators to avoid heap allocations. A logic change to one almost always needs mirroring in the other. Both use the same Bresenham-style incremental bin trackers (`pct_bin`/`abs_bin`) that replace a per-base integer division; the `test_bin_strength_reduction_matches_division` unit test guards them against the direct `floor(i*100/read_len)` / `floor(i/1000)` formulas.
- **`calculate_summary()` is called twice** — once in `main.rs` for console output and again inside `generate_report`. It re-sorts and recomputes each time.
- **Dead code:** `BinnedStatsPercentile::add_read` / `BinnedStatsAbsolute::add_read` are unused; the real accumulation is inlined into `QCAccumulator::add_read*`. Don't assume those methods are the live path.
- Bin dimensions are effectively hard-coded: GC = 101 entries, 2D matrix = `GRID_SIZE` (50) × 50 log-length × linear-quality, percentile = 100 bins, absolute = 101 bins (1 kb × 100 + overflow). The 2D log-length range is `L_MIN=100`..`L_MAX=1_000_000`, quality `0..60`.
- Errors propagate as `Box<dyn std::error::Error>`; `main` prints and `std::process::exit(1)`.
- **The `flate2` direct dependency has the `zlib-ng` feature** purely to feature-unify that backend onto `needletail`'s `flate2`, switching the plain-gzip decode from `miniz_oxide` to zlib-ng (~2–3× faster). It only affects the *serial* plain-`.fastq.gz` fallback — BGZF FASTQ and BAM decode via `noodles` (libdeflate), not `flate2`. Requires `cmake` + a C compiler at build time (zlib-ng builds bundled C via cmake, unlike libdeflater's `cc`-only build). `miniz_oxide` still compiles (needletail keeps `flate2`'s default `rust_backend` feature on) but flate2 selects the zlib-ng C backend when present.
- **The `noodles-bgzf` direct dependency in `Cargo.toml` is not used in code directly** — it exists only to turn on its `libdeflate` feature, which Cargo feature-unifies onto the `noodles-bgzf` instance that the `noodles` meta-crate uses. This switches the BGZF block-inflate backend (`noodles_bgzf::deflate::decode`) from `zlib-rs` to `libdeflater` for faster BAM decompression. It must stay pinned to the same version `noodles` resolves (currently 0.47) or unification breaks and you get two copies. Requires a C compiler at build time (libdeflater builds bundled C).

## Coding guidelines

The `rust-skills` skill (`.agents/skills/rust-skills/`) provides detailed Rust rules; invoke `/rust-skills` when writing or refactoring Rust here.

## Release

`dist-workspace.toml` (cargo-dist) and `recipe/recipe.yaml` (conda / rattler-build, Linux-only) drive packaging via `.github/workflows/release.yml` and `conda.yml`.
