use crate::stats::binned_stats::{BinnedStatsAbsolute, BinnedStatsPercentile};
use serde::Serialize;
use noodles::bam;

const L_MIN: f64 = 100.0;
const L_MAX: f64 = 1_000_000.0; // 1 Mb
const Q_MIN: f64 = 0.0;
const Q_MAX: f64 = 60.0;
pub const GRID_SIZE: usize = 50;

#[derive(Clone, Debug)]
pub struct QCAccumulator {
    pub total_reads: u64,
    pub total_bases: u64,
    pub passed_filters_reads: u64,
    pub passed_filters_bases: u64,

    // Vectors to hold raw read lengths and qualities for precise median & N50 computations.
    pub read_lengths: Vec<u32>,
    pub average_qualities: Vec<f32>,

    pub gc_content_distribution: Vec<u64>, // Size 101: 0% to 100%
    pub length_vs_quality_2d: Vec<Vec<u64>>, // 50x50 grid

    pub quality_by_position_percentile: BinnedStatsPercentile,
    pub base_content_by_position_percentile: BinnedStatsPercentile,

    pub quality_by_position_absolute: BinnedStatsAbsolute,
    pub base_content_by_position_absolute: BinnedStatsAbsolute,
}

#[derive(Debug, Serialize)]
pub struct QCReportSummary {
    pub total_reads: u64,
    pub total_bases: u64,
    pub passed_filters_reads: u64,
    pub passed_filters_bases: u64,
    pub mean_read_length: f64,
    pub median_read_length: u32,
    pub n50_read_length: u32,
    pub min_read_length: u32,
    pub max_read_length: u32,
    pub mean_read_quality: f32,
    pub median_read_quality: f32,
    pub min_read_quality: f32,
    pub max_read_quality: f32,
}

impl QCAccumulator {
    pub fn new() -> Self {
        Self {
            total_reads: 0,
            total_bases: 0,
            passed_filters_reads: 0,
            passed_filters_bases: 0,
            read_lengths: Vec::new(),
            average_qualities: Vec::new(),
            gc_content_distribution: vec![0; 101],
            length_vs_quality_2d: vec![vec![0; GRID_SIZE]; GRID_SIZE],
            quality_by_position_percentile: BinnedStatsPercentile::new(),
            base_content_by_position_percentile: BinnedStatsPercentile::new(),
            quality_by_position_absolute: BinnedStatsAbsolute::new(1000, 100), // 1kb bins up to 100kb+
            base_content_by_position_absolute: BinnedStatsAbsolute::new(1000, 100),
        }
    }

    pub fn add_read(&mut self, seq: &[u8], qual: &[u8], is_fastq_qual: bool, min_len: Option<usize>, min_qual: Option<f32>) {
        let read_len = seq.len();
        self.total_reads += 1;
        self.total_bases += read_len as u64;

        // Calculate average quality
        let mut qual_sum = 0.0;
        for &q_byte in qual.iter() {
            let q = if is_fastq_qual {
                q_byte.saturating_sub(33)
            } else {
                q_byte
            };
            qual_sum += q as f32;
        }
        let avg_qual = if !qual.is_empty() {
            qual_sum / qual.len() as f32
        } else {
            0.0
        };

        // Filter checks
        if let Some(m_len) = min_len {
            if read_len < m_len {
                return;
            }
        }
        if let Some(m_qual) = min_qual {
            if avg_qual < m_qual {
                return;
            }
        }

        self.passed_filters_reads += 1;
        self.passed_filters_bases += read_len as u64;

        // Statistics on passing reads
        self.read_lengths.push(read_len as u32);
        self.average_qualities.push(avg_qual);

        // GC & positional binned stats in a single unified loop
        let mut gc_count = 0;
        let mut valid_bases = 0;
        
        for i in 0..read_len {
            let pct_bin_idx = ((i * 100) / read_len).min(99);
            let abs_bin_idx = (i / 1000).min(100);
            
            let base = seq[i];
            let is_gc = match base {
                b'A' | b'a' => {
                    self.base_content_by_position_percentile.bins[pct_bin_idx].a_count += 1;
                    self.base_content_by_position_absolute.bins[abs_bin_idx].a_count += 1;
                    false
                }
                b'C' | b'c' => {
                    self.base_content_by_position_percentile.bins[pct_bin_idx].c_count += 1;
                    self.base_content_by_position_absolute.bins[abs_bin_idx].c_count += 1;
                    true
                }
                b'G' | b'g' => {
                    self.base_content_by_position_percentile.bins[pct_bin_idx].g_count += 1;
                    self.base_content_by_position_absolute.bins[abs_bin_idx].g_count += 1;
                    true
                }
                b'T' | b't' => {
                    self.base_content_by_position_percentile.bins[pct_bin_idx].t_count += 1;
                    self.base_content_by_position_absolute.bins[abs_bin_idx].t_count += 1;
                    false
                }
                _ => {
                    self.base_content_by_position_percentile.bins[pct_bin_idx].n_count += 1;
                    self.base_content_by_position_absolute.bins[abs_bin_idx].n_count += 1;
                    false
                }
            };
            if is_gc {
                gc_count += 1;
            }
            valid_bases += 1;

            if i < qual.len() {
                let q = if is_fastq_qual {
                    qual[i].saturating_sub(33)
                } else {
                    qual[i]
                };
                let q_val = q as u64;
                self.quality_by_position_percentile.bins[pct_bin_idx].quality_sum += q_val;
                self.quality_by_position_percentile.bins[pct_bin_idx].quality_count += 1;
                self.quality_by_position_absolute.bins[abs_bin_idx].quality_sum += q_val;
                self.quality_by_position_absolute.bins[abs_bin_idx].quality_count += 1;
            }
        }

        let gc_pct = if valid_bases > 0 {
            ((gc_count as f64 / valid_bases as f64) * 100.0).round() as usize
        } else {
            0
        };
        let gc_idx = gc_pct.min(100);
        self.gc_content_distribution[gc_idx] += 1;

        // Length vs Quality 2D Matrix
        let (l_bin, q_bin) = self.get_2d_bin(read_len as u32, avg_qual);
        self.length_vs_quality_2d[l_bin][q_bin] += 1;
    }

    pub fn add_read_bam(&mut self, record: &bam::Record, min_len: Option<usize>, min_qual: Option<f32>) {
        let seq = record.sequence();
        let qual = record.quality_scores();
        let read_len = seq.len();
        
        self.total_reads += 1;
        self.total_bases += read_len as u64;

        // Calculate average quality
        let mut qual_sum = 0.0;
        let mut qual_count = 0;
        for q in qual.iter() {
            qual_sum += q as f32;
            qual_count += 1;
        }
        let avg_qual = if qual_count > 0 {
            qual_sum / qual_count as f32
        } else {
            0.0
        };

        // Filter checks
        if let Some(m_len) = min_len {
            if read_len < m_len {
                return;
            }
        }
        if let Some(m_qual) = min_qual {
            if avg_qual < m_qual {
                return;
            }
        }

        self.passed_filters_reads += 1;
        self.passed_filters_bases += read_len as u64;

        // Statistics on passing reads
        self.read_lengths.push(read_len as u32);
        self.average_qualities.push(avg_qual);

        // GC & positional binned stats in a single unified loop
        let mut gc_count = 0;
        let mut valid_bases = 0;
        
        let mut seq_iter = seq.iter();
        let mut qual_iter = qual.iter();
        
        for i in 0..read_len {
            let pct_bin_idx = ((i * 100) / read_len).min(99);
            let abs_bin_idx = (i / 1000).min(100);
            
            if let Some(base) = seq_iter.next() {
                let is_gc = match base {
                    b'A' | b'a' => {
                        self.base_content_by_position_percentile.bins[pct_bin_idx].a_count += 1;
                        self.base_content_by_position_absolute.bins[abs_bin_idx].a_count += 1;
                        false
                    }
                    b'C' | b'c' => {
                        self.base_content_by_position_percentile.bins[pct_bin_idx].c_count += 1;
                        self.base_content_by_position_absolute.bins[abs_bin_idx].c_count += 1;
                        true
                    }
                    b'G' | b'g' => {
                        self.base_content_by_position_percentile.bins[pct_bin_idx].g_count += 1;
                        self.base_content_by_position_absolute.bins[abs_bin_idx].g_count += 1;
                        true
                    }
                    b'T' | b't' => {
                        self.base_content_by_position_percentile.bins[pct_bin_idx].t_count += 1;
                        self.base_content_by_position_absolute.bins[abs_bin_idx].t_count += 1;
                        false
                    }
                    _ => {
                        self.base_content_by_position_percentile.bins[pct_bin_idx].n_count += 1;
                        self.base_content_by_position_absolute.bins[abs_bin_idx].n_count += 1;
                        false
                    }
                };
                if is_gc {
                    gc_count += 1;
                }
                valid_bases += 1;
            }
            
            if let Some(q) = qual_iter.next() {
                let q_val = q as u64;
                self.quality_by_position_percentile.bins[pct_bin_idx].quality_sum += q_val;
                self.quality_by_position_percentile.bins[pct_bin_idx].quality_count += 1;
                self.quality_by_position_absolute.bins[abs_bin_idx].quality_sum += q_val;
                self.quality_by_position_absolute.bins[abs_bin_idx].quality_count += 1;
            }
        }

        let gc_pct = if valid_bases > 0 {
            ((gc_count as f64 / valid_bases as f64) * 100.0).round() as usize
        } else {
            0
        };
        let gc_idx = gc_pct.min(100);
        self.gc_content_distribution[gc_idx] += 1;

        // Length vs Quality 2D Matrix
        let (l_bin, q_bin) = self.get_2d_bin(read_len as u32, avg_qual);
        self.length_vs_quality_2d[l_bin][q_bin] += 1;
    }

    fn get_2d_bin(&self, length: u32, mean_quality: f32) -> (usize, usize) {
        let l_val = (length as f64).max(L_MIN).min(L_MAX);
        let log_min = L_MIN.ln();
        let log_max = L_MAX.ln();
        let l_pct = (l_val.ln() - log_min) / (log_max - log_min);
        let l_bin = (l_pct * (GRID_SIZE as f64)) as usize;
        let l_bin = l_bin.min(GRID_SIZE - 1);

        let q_val = (mean_quality as f64).max(Q_MIN).min(Q_MAX);
        let q_pct = (q_val - Q_MIN) / (Q_MAX - Q_MIN);
        let q_bin = (q_pct * (GRID_SIZE as f64)) as usize;
        let q_bin = q_bin.min(GRID_SIZE - 1);

        (l_bin, q_bin)
    }

    pub fn merge(&mut self, other: Self) {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.passed_filters_reads += other.passed_filters_reads;
        self.passed_filters_bases += other.passed_filters_bases;

        self.read_lengths.extend(other.read_lengths);
        self.average_qualities.extend(other.average_qualities);

        for i in 0..101 {
            self.gc_content_distribution[i] += other.gc_content_distribution[i];
        }

        for r in 0..GRID_SIZE {
            for c in 0..GRID_SIZE {
                self.length_vs_quality_2d[r][c] += other.length_vs_quality_2d[r][c];
            }
        }

        self.quality_by_position_percentile.merge(&other.quality_by_position_percentile);
        self.base_content_by_position_percentile.merge(&other.base_content_by_position_percentile);
        self.quality_by_position_absolute.merge(&other.quality_by_position_absolute);
        self.base_content_by_position_absolute.merge(&other.base_content_by_position_absolute);
    }

    pub fn calculate_summary(&mut self) -> QCReportSummary {
        let total_reads = self.total_reads;
        let total_bases = self.total_bases;
        let passed_filters_reads = self.passed_filters_reads;
        let passed_filters_bases = self.passed_filters_bases;

        if self.read_lengths.is_empty() {
            return QCReportSummary {
                total_reads,
                total_bases,
                passed_filters_reads,
                passed_filters_bases,
                mean_read_length: 0.0,
                median_read_length: 0,
                n50_read_length: 0,
                min_read_length: 0,
                max_read_length: 0,
                mean_read_quality: 0.0,
                median_read_quality: 0.0,
                min_read_quality: 0.0,
                max_read_quality: 0.0,
            };
        }

        // Mean Read Length
        let sum_len: u64 = self.read_lengths.iter().map(|&x| x as u64).sum();
        let mean_read_length = sum_len as f64 / self.read_lengths.len() as f64;

        // Min/Max Length
        let min_read_length = *self.read_lengths.iter().min().unwrap_or(&0);
        let max_read_length = *self.read_lengths.iter().max().unwrap_or(&0);

        // N50 Read Length
        self.read_lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending
        let half_len = sum_len / 2;
        let mut running_sum = 0u64;
        let mut n50_read_length = self.read_lengths[0];
        for &len in self.read_lengths.iter() {
            running_sum += len as u64;
            if running_sum >= half_len {
                n50_read_length = len;
                break;
            }
        }

        // Median Read Length
        self.read_lengths.sort_unstable(); // Sort ascending
        let mid_l = self.read_lengths.len() / 2;
        let median_read_length = if self.read_lengths.len() % 2 == 0 {
            ((self.read_lengths[mid_l - 1] as u64 + self.read_lengths[mid_l] as u64) / 2) as u32
        } else {
            self.read_lengths[mid_l]
        };

        // Mean Quality
        let sum_q: f32 = self.average_qualities.iter().sum();
        let mean_read_quality = sum_q / self.average_qualities.len() as f32;

        // Min/Max Quality
        let min_read_quality = self.average_qualities.iter().copied().fold(f32::INFINITY, f32::min);
        let max_read_quality = self.average_qualities.iter().copied().fold(f32::NEG_INFINITY, f32::max);

        // Median Quality
        self.average_qualities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid_q = self.average_qualities.len() / 2;
        let median_read_quality = if self.average_qualities.len() % 2 == 0 {
            (self.average_qualities[mid_q - 1] + self.average_qualities[mid_q]) / 2.0
        } else {
            self.average_qualities[mid_q]
        };

        QCReportSummary {
            total_reads,
            total_bases,
            passed_filters_reads,
            passed_filters_bases,
            mean_read_length,
            median_read_length,
            n50_read_length,
            min_read_length,
            max_read_length,
            mean_read_quality,
            median_read_quality,
            min_read_quality,
            max_read_quality,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n50_and_median() {
        let mut acc = QCAccumulator::new();
        // Add 5 reads of lengths: 10, 20, 30, 40, 100
        // Total bases = 200. N50 is the length where cumulative sum >= 100
        // Descending order of lengths: 100, 40, 30, 20, 10
        // Cumulative sums: 100 (which is >= 100), so N50 should be 100!
        // Median of 5 reads is 30.
        acc.add_read(b"A".repeat(10).as_slice(), vec![30; 10].as_slice(), false, None, None);
        acc.add_read(b"A".repeat(20).as_slice(), vec![30; 20].as_slice(), false, None, None);
        acc.add_read(b"A".repeat(30).as_slice(), vec![30; 30].as_slice(), false, None, None);
        acc.add_read(b"A".repeat(40).as_slice(), vec![30; 40].as_slice(), false, None, None);
        acc.add_read(b"A".repeat(100).as_slice(), vec![30; 100].as_slice(), false, None, None);

        let summary = acc.calculate_summary();
        assert_eq!(summary.total_reads, 5);
        assert_eq!(summary.total_bases, 200);
        assert_eq!(summary.n50_read_length, 100);
        assert_eq!(summary.median_read_length, 30);
        assert_eq!(summary.mean_read_length, 40.0);
    }

    #[test]
    fn test_filters() {
        let mut acc = QCAccumulator::new();
        // Add a read with length 10, mean quality 20 (raw Phred)
        acc.add_read(b"A".repeat(10).as_slice(), vec![20; 10].as_slice(), false, Some(15), None); // fails length filter
        acc.add_read(b"A".repeat(20).as_slice(), vec![10; 20].as_slice(), false, None, Some(15.0)); // fails quality filter
        acc.add_read(b"A".repeat(20).as_slice(), vec![20; 20].as_slice(), false, Some(15), Some(15.0)); // passes both

        let summary = acc.calculate_summary();
        assert_eq!(summary.total_reads, 3);
        assert_eq!(summary.passed_filters_reads, 1);
        assert_eq!(summary.passed_filters_bases, 20);
    }
}
