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
    
    // Percentile binned data (0% to 100% of read length)
    pub quality_by_position_percentile: Vec<f32>,
    pub base_content_by_position_percentile: PercentileBaseContent,
    
    // Absolute binned data (e.g. 1kb intervals)
    pub quality_by_position_absolute: Vec<f32>,
    pub base_content_by_position_absolute: AbsoluteBaseContent,
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

pub fn generate_report(mut accumulator: QCAccumulator, output_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let summary = accumulator.calculate_summary();
    
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
