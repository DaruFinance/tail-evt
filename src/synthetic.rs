//! Deterministic synthetic-returns generator used by `--synthetic`.
//!
//! The four series are built from Student-t(ν=5) draws, scaled to a daily-return
//! magnitude. By Pickands–Balkema–de Haan applied to the t-distribution, the
//! upper tail of t(ν) is asymptotically GPD with shape ξ = 1/ν, so a peaks-
//! over-threshold fit on this data should recover ξ close to 0.2.
//!
//! Cross-asset structure:
//!   A_lead  = F + ε₁         shared factor F creates non-trivial χ
//!   A_lag   = 0.7·F + ε₂
//!   B_indep1 = ε₃            independent
//!   B_indep2 = ε₄            independent
//!
//! All ε's and F are independent t(5)·σ draws; deterministic via SEED.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, StudentT};

pub const ASSETS: [&str; 4] = ["A_lead", "A_lag", "B_indep1", "B_indep2"];
pub const SEED: u64 = 2026_04_25;

/// Degrees of freedom for the Student-t marginal noise. ν = 5 ⇒ tail shape
/// parameter ξ ≈ 1/ν = 0.20 in the GPD limit.
pub const T_DOF: f64 = 5.0;

/// Per-component scale: a t(5) draw scaled by SCALE has the rough magnitude
/// of a daily equity return.
pub const SCALE: f64 = 0.004;

pub fn generate(n: usize) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let t = StudentT::new(T_DOF).expect("valid degrees of freedom");

    let mut a_lead = Vec::with_capacity(n);
    let mut a_lag = Vec::with_capacity(n);
    let mut b_1 = Vec::with_capacity(n);
    let mut b_2 = Vec::with_capacity(n);

    for _ in 0..n {
        let f = t.sample(&mut rng) * SCALE;
        let e1 = t.sample(&mut rng) * SCALE;
        let e2 = t.sample(&mut rng) * SCALE;
        let e3 = t.sample(&mut rng) * SCALE;
        let e4 = t.sample(&mut rng) * SCALE;

        a_lead.push(f + e1);
        a_lag.push(0.7 * f + e2);
        b_1.push(e3);
        b_2.push(e4);
    }
    // (avoid an unused-warning on rand::Rng — kept for future use)
    let _ = rng.gen::<u8>();

    vec![a_lead, a_lag, b_1, b_2]
}

pub fn labels() -> Vec<String> {
    ASSETS.iter().map(|s| s.to_string()).collect()
}
