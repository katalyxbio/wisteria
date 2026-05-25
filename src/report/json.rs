use crate::stats::{QCAccumulator, QCReportSummary};
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Serialize)]
pub struct UnifiedQCReport {
    pub summary: QCReportSummary,
    pub gc_content_distribution: Vec<u64>,
    pub length_vs_quality_2d: Vec<Vec<u64>>,
    
    // Explicit 1D Distributions for plotting
    pub sequence_length_distribution: Histogram1D<u32>,
    pub per_sequence_quality_distribution: Histogram1D<f32>,
    
    // Percentile binned data (0% to 100% of read length)
    pub quality_by_position_percentile: Vec<f32>,
    pub base_content_by_position_percentile: PercentileBaseContent,
    
    // Absolute binned data (e.g. 1kb intervals)
    pub quality_by_position_absolute: Vec<f32>,
    pub base_content_by_position_absolute: AbsoluteBaseContent,
}

#[derive(Serialize)]
pub struct Histogram1D<T> {
    pub bin_edges: Vec<T>,
    pub counts: Vec<u64>,
}

#[derive(Serialize)]
pub struct PercentileBaseContent {
    pub a: Vec<u64>,
    pub c: Vec<u64>,
    pub g: Vec<u64>,
    pub t: Vec<u64>,
    pub n: Vec<u64>,
}

#[derive(Serialize)]
pub struct AbsoluteBaseContent {
    pub a: Vec<u64>,
    pub c: Vec<u64>,
    pub g: Vec<u64>,
    pub t: Vec<u64>,
    pub n: Vec<u64>,
}

fn bin_read_lengths(lengths: &[u32]) -> Histogram1D<u32> {
    if lengths.is_empty() {
        return Histogram1D { bin_edges: vec![0; 51], counts: vec![0; 50] };
    }
    
    let min_l = *lengths.iter().min().unwrap() as f64;
    let max_l = *lengths.iter().max().unwrap() as f64;
    
    let bin_count = 50;
    let mut counts = vec![0u64; bin_count];
    let mut bin_edges = Vec::with_capacity(bin_count + 1);
    
    let step = if max_l > min_l {
        (max_l - min_l) / bin_count as f64
    } else {
        1.0
    };
    
    for i in 0..=bin_count {
        bin_edges.push((min_l + i as f64 * step).round() as u32);
    }
    
    for &len in lengths {
        let val = len as f64;
        let mut bin = if step > 0.0 {
            ((val - min_l) / step) as usize
        } else {
            0
        };
        if bin >= bin_count {
            bin = bin_count - 1;
        }
        counts[bin] += 1;
    }
    
    Histogram1D { bin_edges, counts }
}

fn bin_average_qualities(quals: &[f32]) -> Histogram1D<f32> {
    let bin_count = 60;
    let mut counts = vec![0u64; bin_count];
    let mut bin_edges = Vec::with_capacity(bin_count + 1);
    
    for i in 0..=bin_count {
        bin_edges.push(i as f32);
    }
    
    for &q in quals {
        let mut bin = q.floor() as usize;
        if bin >= bin_count {
            bin = bin_count - 1;
        }
        counts[bin] += 1;
    }
    
    Histogram1D { bin_edges, counts }
}

pub fn generate_report(mut accumulator: QCAccumulator, output_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let summary = accumulator.calculate_summary();
    
    let sequence_length_distribution = bin_read_lengths(&accumulator.read_lengths);
    let per_sequence_quality_distribution = bin_average_qualities(&accumulator.average_qualities);
    
    let quality_by_position_percentile: Vec<f32> = accumulator.quality_by_position_percentile.bins.iter()
        .map(|bin| {
            if bin.quality_count > 0 {
                bin.quality_sum as f32 / bin.quality_count as f32
            } else {
                0.0
            }
        })
        .collect();
        
    let mut a_pct = Vec::with_capacity(100);
    let mut c_pct = Vec::with_capacity(100);
    let mut g_pct = Vec::with_capacity(100);
    let mut t_pct = Vec::with_capacity(100);
    let mut n_pct = Vec::with_capacity(100);
    
    for bin in accumulator.base_content_by_position_percentile.bins.iter() {
        a_pct.push(bin.a_count);
        c_pct.push(bin.c_count);
        g_pct.push(bin.g_count);
        t_pct.push(bin.t_count);
        n_pct.push(bin.n_count);
    }
    
    let base_content_by_position_percentile = PercentileBaseContent {
        a: a_pct,
        c: c_pct,
        g: g_pct,
        t: t_pct,
        n: n_pct,
    };
    
    let quality_by_position_absolute: Vec<f32> = accumulator.quality_by_position_absolute.bins.iter()
        .map(|bin| {
            if bin.quality_count > 0 {
                bin.quality_sum as f32 / bin.quality_count as f32
            } else {
                0.0
            }
        })
        .collect();
        
    let mut a_abs = Vec::new();
    let mut c_abs = Vec::new();
    let mut g_abs = Vec::new();
    let mut t_abs = Vec::new();
    let mut n_abs = Vec::new();
    
    for bin in accumulator.base_content_by_position_absolute.bins.iter() {
        a_abs.push(bin.a_count);
        c_abs.push(bin.c_count);
        g_abs.push(bin.g_count);
        t_abs.push(bin.t_count);
        n_abs.push(bin.n_count);
    }
    
    let base_content_by_position_absolute = AbsoluteBaseContent {
        a: a_abs,
        c: c_abs,
        g: g_abs,
        t: t_abs,
        n: n_abs,
    };
    
    let report = UnifiedQCReport {
        summary,
        gc_content_distribution: accumulator.gc_content_distribution,
        length_vs_quality_2d: accumulator.length_vs_quality_2d,
        sequence_length_distribution,
        per_sequence_quality_distribution,
        quality_by_position_percentile,
        base_content_by_position_percentile,
        quality_by_position_absolute,
        base_content_by_position_absolute,
    };
    
    let serialized = serde_json::to_string_pretty(&report)?;
    let mut file = File::create(output_path)?;
    file.write_all(serialized.as_bytes())?;
    
    Ok(())
}
