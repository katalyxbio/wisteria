//! PNG plot rendering for the `--plots` flag.
//!
//! Renders the [`UnifiedQCReport`] into a set of standalone PNG charts using
//! `plotters` (pure-Rust: bitmap backend + `ab_glyph` fonts, no C deps). One
//! file per chart is written into the target directory. Kept independent of the
//! JSON path so the two outputs don't interfere.

use super::json::UnifiedQCReport;
use plotters::prelude::*;
use plotters::style::{register_font, FontStyle};
use std::path::Path;
use std::sync::Once;

/// Supersample factor. The pure-Rust `ab_glyph` text layout rounds each glyph's
/// pen position to a whole pixel (no subpixel positioning/kerning), so at small
/// font sizes advance-width rounding shows up as uneven letter spacing and
/// slightly off-center captions. Rendering at a higher pixel density with
/// proportionally larger fonts makes each glyph advance span more pixels, so the
/// ±0.5 px rounding is a smaller fraction of the advance and the text smooths
/// out (and the PNGs come out crisper). All layout sizes below scale with this.
const SCALE: u32 = 2;

/// Canvas size for every chart (px).
const W: u32 = 960 * SCALE;
const H: u32 = 600 * SCALE;

// Layout sizes, all scaled so proportions stay identical as SCALE changes.
const CAPTION_PT: i32 = 26 * SCALE as i32;
const LABEL_PT: i32 = 15 * SCALE as i32;
const LEGEND_PT: i32 = 16 * SCALE as i32;
const MARGIN: i32 = 16 * SCALE as i32;
const X_AXIS: i32 = 52 * SCALE as i32;
const Y_AXIS: i32 = 72 * SCALE as i32;
const STROKE: u32 = 2 * SCALE;

/// Bundled font. The pure-Rust `ab_glyph` backend ships no default font, so we
/// vendor one (Roboto, Apache-2.0) and register it under the family names the
/// charts request. Without this, rendering fails with `FontUnavailable`.
static FONT_BYTES: &[u8] = include_bytes!("../../resources/Roboto-Regular.ttf");
static FONT_INIT: Once = Once::new();

fn ensure_font() {
    FONT_INIT.call_once(|| {
        // Register for every family name used in captions/labels below.
        for family in ["sans-serif", "serif", "monospace"] {
            let _ = register_font(family, FontStyle::Normal, FONT_BYTES);
        }
    });
}

type PlotResult = Result<(), Box<dyn std::error::Error>>;

/// Render the full set of QC charts as PNGs into `dir` (created if absent).
///
/// Returns the list of files written so the caller can report them.
pub fn generate_plots(report: &UnifiedQCReport, dir: &Path) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    ensure_font();
    std::fs::create_dir_all(dir)?;
    let mut written = Vec::new();

    let mut render = |name: &str, f: &dyn Fn(&Path) -> PlotResult| -> Result<(), Box<dyn std::error::Error>> {
        let path = dir.join(name);
        f(&path)?;
        written.push(name.to_string());
        Ok(())
    };

    render("read_length_distribution.png", &|p| plot_length_hist(report, p))?;
    render("per_sequence_quality_distribution.png", &|p| plot_quality_hist(report, p))?;
    render("gc_content_distribution.png", &|p| plot_gc(report, p))?;
    render("quality_by_position_percentile.png", &|p| {
        plot_quality_by_position(
            &report.quality_by_position_percentile,
            "Mean Quality by Position (relative %)",
            "Read position (% of length)",
            1.0,
            p,
        )
    })?;
    render("quality_by_position_absolute.png", &|p| {
        plot_quality_by_position(
            &report.quality_by_position_absolute,
            "Mean Quality by Position (absolute)",
            "Read position (kb bin, last = 100kb+)",
            1.0,
            p,
        )
    })?;
    render("base_content_by_position.png", &|p| plot_base_content(report, p))?;
    render("length_vs_quality_heatmap.png", &|p| plot_heatmap(report, p))?;

    Ok(written)
}

/// Draw filled bars from `(x_left, x_right, height)` triples.
fn draw_bars(
    path: &Path,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
    bars: &[(f64, f64, f64)],
    color: RGBColor,
) -> PlotResult {
    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_min = bars.first().map(|b| b.0).unwrap_or(0.0);
    let x_max = bars.last().map(|b| b.1).unwrap_or(1.0).max(x_min + 1.0);
    let y_max = bars.iter().map(|b| b.2).fold(0.0_f64, f64::max).max(1.0) * 1.05;

    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", CAPTION_PT))
        .margin(MARGIN)
        .x_label_area_size(X_AXIS)
        .y_label_area_size(Y_AXIS)
        .build_cartesian_2d(x_min..x_max, 0.0..y_max)?;

    chart
        .configure_mesh()
        .x_desc(x_desc)
        .y_desc(y_desc)
        .label_style(("sans-serif", LABEL_PT))
        .draw()?;

    chart.draw_series(bars.iter().map(|&(x0, x1, h)| {
        Rectangle::new([(x0, 0.0), (x1, h)], color.mix(0.85).filled())
    }))?;

    root.present()?;
    Ok(())
}

fn plot_length_hist(report: &UnifiedQCReport, path: &Path) -> PlotResult {
    let hist = &report.sequence_length_distribution;
    let edges = &hist.bin_edges;
    let bars: Vec<(f64, f64, f64)> = hist
        .counts
        .iter()
        .enumerate()
        .map(|(i, &c)| (edges[i] as f64, edges[i + 1] as f64, c as f64))
        .collect();
    draw_bars(
        path,
        "Read Length Distribution",
        "Read length (bp)",
        "Read count",
        &bars,
        BLUE,
    )
}

fn plot_quality_hist(report: &UnifiedQCReport, path: &Path) -> PlotResult {
    let hist = &report.per_sequence_quality_distribution;
    let edges = &hist.bin_edges;
    let bars: Vec<(f64, f64, f64)> = hist
        .counts
        .iter()
        .enumerate()
        .map(|(i, &c)| (edges[i] as f64, edges[i + 1] as f64, c as f64))
        .collect();
    draw_bars(
        path,
        "Per-Sequence Quality Distribution",
        "Mean Phred quality (Q)",
        "Read count",
        &bars,
        RGBColor(0x2e, 0x8b, 0x57),
    )
}

fn plot_gc(report: &UnifiedQCReport, path: &Path) -> PlotResult {
    let bars: Vec<(f64, f64, f64)> = report
        .gc_content_distribution
        .iter()
        .enumerate()
        .map(|(i, &c)| (i as f64, i as f64 + 1.0, c as f64))
        .collect();
    draw_bars(
        path,
        "GC Content Distribution",
        "GC content (%)",
        "Read count",
        &bars,
        RGBColor(0xb2, 0x22, 0x22),
    )
}

/// Line chart of mean quality across position bins. `x_step` scales the bin
/// index onto the x-axis (1.0 keeps raw bin indices).
fn plot_quality_by_position(
    values: &[f32],
    caption: &str,
    x_desc: &str,
    x_step: f64,
    path: &Path,
) -> PlotResult {
    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_max = (values.len() as f64 * x_step).max(1.0);
    let y_max = values.iter().cloned().fold(0.0_f32, f32::max).max(1.0) as f64 * 1.1;

    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", CAPTION_PT))
        .margin(MARGIN)
        .x_label_area_size(X_AXIS)
        .y_label_area_size(Y_AXIS)
        .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;

    chart
        .configure_mesh()
        .x_desc(x_desc)
        .y_desc("Mean Phred quality (Q)")
        .label_style(("sans-serif", LABEL_PT))
        .draw()?;

    chart.draw_series(LineSeries::new(
        values
            .iter()
            .enumerate()
            .map(|(i, &q)| (i as f64 * x_step, q as f64)),
        RGBColor(0x1f, 0x77, 0xb4).stroke_width(STROKE),
    ))?;

    root.present()?;
    Ok(())
}

fn plot_base_content(report: &UnifiedQCReport, path: &Path) -> PlotResult {
    let bc = &report.base_content_by_position_percentile;
    let n = bc.a.len();

    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Base Content by Position (relative %)", ("sans-serif", CAPTION_PT))
        .margin(MARGIN)
        .x_label_area_size(X_AXIS)
        .y_label_area_size(Y_AXIS)
        .build_cartesian_2d(0.0..(n as f64).max(1.0), 0.0..100.0)?;

    chart
        .configure_mesh()
        .x_desc("Read position (% of length)")
        .y_desc("Base fraction (%)")
        .label_style(("sans-serif", LABEL_PT))
        .draw()?;

    // Each channel as its own line; y is the per-position percentage of that base.
    let channels: [(&str, &Vec<u64>, RGBColor); 5] = [
        ("A", &bc.a, RGBColor(0x2c, 0xa0, 0x2c)), // green
        ("C", &bc.c, RGBColor(0x1f, 0x77, 0xb4)), // blue
        ("G", &bc.g, RGBColor(0xff, 0x7f, 0x0e)), // orange
        ("T", &bc.t, RGBColor(0xd6, 0x27, 0x28)), // red
        ("N", &bc.n, RGBColor(0x7f, 0x7f, 0x7f)), // gray
    ];

    for (label, counts, color) in channels {
        let series = (0..n).map(|i| {
            let total = bc.a[i] + bc.c[i] + bc.g[i] + bc.t[i] + bc.n[i];
            let pct = if total > 0 {
                counts[i] as f64 / total as f64 * 100.0
            } else {
                0.0
            };
            (i as f64, pct)
        });
        chart
            .draw_series(LineSeries::new(series, color.stroke_width(STROKE)))?
            .label(label)
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], color.stroke_width(STROKE + 1)));
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK.mix(0.3))
        .label_font(("sans-serif", LEGEND_PT))
        .draw()?;

    root.present()?;
    Ok(())
}

/// Log-count color ramp (light → dark blue) for the density heatmap.
fn heat_color(t: f64) -> RGBColor {
    // t in [0,1]; blend near-white to deep blue-purple.
    let t = t.clamp(0.0, 1.0);
    let r = (245.0 - t * (245.0 - 40.0)) as u8;
    let g = (245.0 - t * (245.0 - 20.0)) as u8;
    let b = (250.0 - t * (250.0 - 130.0)) as u8;
    RGBColor(r, g, b)
}

fn plot_heatmap(report: &UnifiedQCReport, path: &Path) -> PlotResult {
    let grid = &report.length_vs_quality_2d;
    let nx = grid.len();
    let ny = grid.first().map(|r| r.len()).unwrap_or(0);
    if nx == 0 || ny == 0 {
        // Nothing to draw; emit an empty canvas so the file still exists.
        let root = BitMapBackend::new(path, (W, H)).into_drawing_area();
        root.fill(&WHITE)?;
        root.present()?;
        return Ok(());
    }

    // Log-scale the counts so dense clusters don't wash everything else out.
    let max_count = grid
        .iter()
        .flat_map(|row| row.iter())
        .cloned()
        .max()
        .unwrap_or(0);
    let log_max = ((max_count as f64) + 1.0).ln().max(1.0);

    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Read Length vs. Quality Density", ("sans-serif", CAPTION_PT))
        .margin(MARGIN)
        .x_label_area_size(X_AXIS)
        .y_label_area_size(Y_AXIS)
        .build_cartesian_2d(0.0..nx as f64, 0.0..ny as f64)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("Read length bin (log 100 bp – 1 Mb)")
        .y_desc("Mean quality bin (Q0 – Q60)")
        .label_style(("sans-serif", LABEL_PT))
        .draw()?;

    let mut cells = Vec::with_capacity(nx * ny);
    for (xi, row) in grid.iter().enumerate() {
        for (yi, &count) in row.iter().enumerate() {
            let t = if count > 0 {
                ((count as f64) + 1.0).ln() / log_max
            } else {
                0.0
            };
            cells.push(Rectangle::new(
                [(xi as f64, yi as f64), (xi as f64 + 1.0, yi as f64 + 1.0)],
                heat_color(t).filled(),
            ));
        }
    }
    chart.draw_series(cells)?;

    root.present()?;
    Ok(())
}
