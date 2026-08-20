pub mod bam;
pub mod fastq;

pub use bam::{benchmark_decompress, process_bam};
pub use fastq::process_fastq;
