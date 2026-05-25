use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct PositionStats {
    pub quality_sum: u64,
    pub quality_count: u64,
    pub a_count: u64,
    pub c_count: u64,
    pub g_count: u64,
    pub t_count: u64,
    pub n_count: u64,
}

impl PositionStats {
    pub fn new() -> Self {
        Self {
            quality_sum: 0,
            quality_count: 0,
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            n_count: 0,
        }
    }

    pub fn merge(&mut self, other: &Self) {
        self.quality_sum += other.quality_sum;
        self.quality_count += other.quality_count;
        self.a_count += other.a_count;
        self.c_count += other.c_count;
        self.g_count += other.g_count;
        self.t_count += other.t_count;
        self.n_count += other.n_count;
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct BinnedStatsPercentile {
    pub bins: Vec<PositionStats>,
}

impl BinnedStatsPercentile {
    pub fn new() -> Self {
        Self {
            bins: vec![PositionStats::new(); 100],
        }
    }

    pub fn merge(&mut self, other: &Self) {
        for (i, bin) in self.bins.iter_mut().enumerate() {
            bin.merge(&other.bins[i]);
        }
    }

    pub fn add_read(&mut self, seq: &[u8], qual: &[u8], is_fastq_qual: bool) {
        let len = seq.len();
        if len == 0 {
            return;
        }

        for i in 0..len {
            let bin_idx = (i * 100) / len;
            let bin_idx = bin_idx.min(99);
            let bin = &mut self.bins[bin_idx];

            let base = seq[i];
            match base {
                b'A' | b'a' => bin.a_count += 1,
                b'C' | b'c' => bin.c_count += 1,
                b'G' | b'g' => bin.g_count += 1,
                b'T' | b't' => bin.t_count += 1,
                _ => bin.n_count += 1,
            }

            if i < qual.len() {
                let q = if is_fastq_qual {
                    qual[i].saturating_sub(33)
                } else {
                    qual[i]
                };
                bin.quality_sum += q as u64;
                bin.quality_count += 1;
            }
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct BinnedStatsAbsolute {
    pub bins: Vec<PositionStats>,
    pub bin_size: usize,
    pub max_bins: usize,
}

impl BinnedStatsAbsolute {
    pub fn new(bin_size: usize, max_bins: usize) -> Self {
        Self {
            bins: vec![PositionStats::new(); max_bins + 1],
            bin_size,
            max_bins,
        }
    }

    pub fn merge(&mut self, other: &Self) {
        for (i, bin) in self.bins.iter_mut().enumerate() {
            if i < other.bins.len() {
                bin.merge(&other.bins[i]);
            }
        }
    }

    pub fn add_read(&mut self, seq: &[u8], qual: &[u8], is_fastq_qual: bool) {
        let len = seq.len();
        if len == 0 {
            return;
        }

        for i in 0..len {
            let bin_idx = i / self.bin_size;
            let bin_idx = bin_idx.min(self.max_bins);
            let bin = &mut self.bins[bin_idx];

            let base = seq[i];
            match base {
                b'A' | b'a' => bin.a_count += 1,
                b'C' | b'c' => bin.c_count += 1,
                b'G' | b'g' => bin.g_count += 1,
                b'T' | b't' => bin.t_count += 1,
                _ => bin.n_count += 1,
            }

            if i < qual.len() {
                let q = if is_fastq_qual {
                    qual[i].saturating_sub(33)
                } else {
                    qual[i]
                };
                bin.quality_sum += q as u64;
                bin.quality_count += 1;
            }
        }
    }
}
