//! tail-evt: peaks-over-threshold GPD fits + pairwise χ tail coupling.
//!
//! Two operating modes:
//!
//!   * `--synthetic`        Generate a deterministic 4-asset synthetic universe
//!                          inline (see `synthetic.rs`) and run the full
//!                          pipeline against it. No external CSV required.
//!
//!   * `--input <CSV>`      Read returns from a CSV file (header row of asset
//!                          names, numeric rows of contemporaneous returns).
//!
//! Output is a single JSON document with per-asset GPD fits (ξ, β, threshold,
//! exceedance count) and a pairwise empirical χ matrix.

mod synthetic;

use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use std::{fs, path::PathBuf};

#[derive(Parser, Debug)]
#[command(about = "Peaks-over-threshold GPD fits + pairwise χ tail coupling")]
struct Args {
    /// CSV file of returns (header = asset names, numeric rows = timesteps).
    /// Required unless --synthetic is passed.
    #[arg(long)]
    input: Option<PathBuf>,

    /// Run the deterministic synthetic-returns demo instead of reading a CSV.
    #[arg(long, default_value_t = false)]
    synthetic: bool,

    /// Number of timesteps in the synthetic demo (ignored when --input is given).
    #[arg(long, default_value_t = 6000)]
    synthetic_n: usize,

    /// Upper-tail quantile used as the GPD threshold (e.g. 0.95).
    #[arg(long, default_value_t = 0.95)]
    threshold_quantile: f64,

    /// Output JSON path.
    #[arg(long, default_value = "tail.json")]
    output: PathBuf,
}

#[derive(Serialize)]
struct AssetFit {
    asset: String,
    u: f64,
    xi: f64,
    beta: f64,
    n_exc: usize,
}

#[derive(Serialize)]
struct Output {
    fits: Vec<AssetFit>,
    labels: Vec<String>,
    chi_matrix: Vec<Vec<f64>>,
    threshold_quantile: f64,
    n_obs: usize,
    source: &'static str,
}

fn read_csv(path: &PathBuf) -> (Vec<String>, Vec<Vec<f64>>) {
    let txt = fs::read_to_string(path).expect("read CSV");
    let mut lines = txt.lines();
    let header: Vec<String> = lines
        .next()
        .expect("empty CSV")
        .split(',')
        .map(str::to_string)
        .collect();
    let mut cols: Vec<Vec<f64>> = (0..header.len()).map(|_| Vec::new()).collect();
    for line in lines {
        for (j, cell) in line.split(',').enumerate() {
            cols[j].push(cell.trim().parse().unwrap_or(f64::NAN));
        }
    }
    (header, cols)
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let idx = ((n as f64 - 1.0) * q).round() as usize;
    sorted[idx.min(n - 1)]
}

/// Probability-weighted-moments estimator of GPD(ξ, β) (Hosking & Wallis 1987).
///
/// Sample probability-weighted moments, with x_(i) the i-th order statistic
/// (0-based, ascending), are
///   β̂₀ = (1/n) Σ x_(i)
///   β̂₁ = (1/n) Σ (i/(n−1)) · x_(i)        // weight is the empirical F̂(x_(i))
///
/// For the GPD with shape ξ and scale σ these have the closed forms
///   β₀ = σ/(1−ξ)
///   β₁ = σ/(1−ξ) − σ/(2(2−ξ))
///
/// Inverting the system (and noting the sign flip relative to the form quoted
/// for the (1−F̂)-weighted PWM) gives
///   ξ = 2 + β̂₀ / (β̂₀ − 2 β̂₁)
///   σ = β̂₀ · (1 − ξ)
///
/// ξ is clamped to [−0.9, 0.9] and σ is floored at 1e-8 for numerical safety
/// on small or near-degenerate samples. A direct sanity check: for an
/// exponential population (ξ = 0) the closed form gives β₁ = 3β/4, hence
/// β̂₀ / (β̂₀ − 2 β̂₁) → 1/(1 − 3/2) = −2, and ξ̂ → 2 + (−2) = 0.
fn fit_gpd(excesses: &[f64]) -> (f64, f64) {
    if excesses.len() < 30 {
        return (
            0.0,
            excesses.iter().sum::<f64>() / excesses.len().max(1) as f64,
        );
    }
    let mut s = excesses.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len() as f64;
    let b0 = s.iter().sum::<f64>() / n;
    let b1 = s
        .iter()
        .enumerate()
        .map(|(i, &x)| (i as f64 / (n - 1.0)) * x)
        .sum::<f64>()
        / n;
    let denom = b0 - 2.0 * b1;
    if denom.abs() < 1e-12 {
        return (0.0, b0.max(1e-8));
    }
    let xi = (2.0 + b0 / denom).clamp(-0.9, 0.9);
    let beta = (b0 * (1.0 - xi)).max(1e-8);
    (xi, beta)
}

/// Empirical χ(u) = P(Y > u_Y | X > u_X), with u_X / u_Y matched-quantile thresholds.
fn chi(x: &[f64], y: &[f64], ux: f64, uy: f64) -> f64 {
    let mut num = 0usize;
    let mut den = 0usize;
    for (a, b) in x.iter().zip(y) {
        if a.is_finite() && *a > ux {
            den += 1;
            if b.is_finite() && *b > uy {
                num += 1;
            }
        }
    }
    if den == 0 { 0.0 } else { num as f64 / den as f64 }
}

fn main() {
    let args = Args::parse();

    let (labels, cols, source): (Vec<String>, Vec<Vec<f64>>, &'static str) = match (&args.input, args.synthetic) {
        (Some(path), _) => {
            let (h, c) = read_csv(path);
            (h, c, "csv")
        }
        (None, true) => (synthetic::labels(), synthetic::generate(args.synthetic_n), "synthetic"),
        (None, false) => {
            eprintln!("error: pass --input <CSV> or --synthetic.");
            std::process::exit(2);
        }
    };

    let n_obs = cols.first().map(|c| c.len()).unwrap_or(0);

    // per-asset threshold + GPD fit, in parallel
    let per_asset: Vec<(String, f64, (f64, f64), Vec<f64>)> = labels
        .par_iter()
        .enumerate()
        .map(|(j, name)| {
            let col = &cols[j];
            let mut s: Vec<f64> = col.iter().copied().filter(|v| v.is_finite()).collect();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let u = quantile(&s, args.threshold_quantile);
            let exc: Vec<f64> = col
                .iter()
                .copied()
                .filter(|v| v.is_finite() && *v > u)
                .map(|v| v - u)
                .collect();
            let fit = fit_gpd(&exc);
            (name.clone(), u, fit, exc)
        })
        .collect();

    let fits: Vec<AssetFit> = per_asset
        .iter()
        .map(|(n, u, (xi, b), e)| AssetFit {
            asset: n.clone(),
            u: *u,
            xi: *xi,
            beta: *b,
            n_exc: e.len(),
        })
        .collect();

    // pairwise χ matrix
    let k = labels.len();
    let mut chi_m = vec![vec![0.0f64; k]; k];
    for i in 0..k {
        for j in 0..k {
            if i == j {
                chi_m[i][j] = 1.0;
                continue;
            }
            let ui = per_asset[i].1;
            let uj = per_asset[j].1;
            chi_m[i][j] = chi(&cols[i], &cols[j], ui, uj);
        }
    }

    let out = Output {
        fits,
        labels,
        chi_matrix: chi_m,
        threshold_quantile: args.threshold_quantile,
        n_obs,
        source,
    };
    fs::write(&args.output, serde_json::to_string_pretty(&out).unwrap()).unwrap();
    eprintln!("wrote {}", args.output.display());
}
