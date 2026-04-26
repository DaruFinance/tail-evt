# tail-evt

**Extreme-value analysis of cross-asset tail dependence.**

Peaks-over-threshold with Generalised Pareto fits on joint asset exceedances.
A compact Rust core handles per-asset GPD fitting and pairwise empirical χ
estimation in parallel (Rayon). The repository ships a deterministic synthetic
demo so `git clone && cargo run --release -- --synthetic` produces real
output without any external data.

## Reproduce

```bash
git clone https://github.com/DaruFinance/tail-evt
cd tail-evt
cargo run --release -- --synthetic --output tail_synthetic.json

# (optional) render demo PNGs from the JSON
python scripts/render_demo.py tail_synthetic.json figures/
```

## Problem statement

For a return series `Xₜ` and threshold `u`, the excesses `Xₜ − u | Xₜ > u`
converge (Pickands–Balkema–de Haan) to the Generalised Pareto distribution

```
G_{ξ, β}(y) = 1 − (1 + ξ y / β)^{−1/ξ},   y ≥ 0
```

with shape `ξ` and scale `β`. For a pair of assets `(X, Y)`, tail coupling
is summarised by the χ-coefficient

```
χ = lim_{u → F^{−1}(1)} P(Y > u | X > u)
```

estimated via a consistent empirical plug-in. This repository estimates
`(ξ, β)` per-asset via the probability-weighted-moments estimator
(Hosking & Wallis 1987) and `χ(u)` per-pair at the same threshold quantile,
returning the joint tail-dependence matrix.

## Usage

```bash
# build the release binary
cargo build --release

# Mode A: deterministic synthetic universe (4 assets, no external data)
./target/release/evt-cli --synthetic \
    --synthetic-n 6000 \
    --threshold-quantile 0.95 \
    --output tail_synthetic.json

# Mode B: real CSV (header = asset names, numeric rows = contemporaneous returns)
./target/release/evt-cli \
    --input returns.csv \
    --threshold-quantile 0.95 \
    --output tail.json
```

Either mode emits a JSON document with:

- `fits` — per-asset `(asset, u, ξ, β, n_exc)`
- `labels` — asset order
- `chi_matrix` — pairwise empirical χ(u) at the chosen threshold
- `threshold_quantile`, `n_obs`, `source` (`"synthetic"` or `"csv"`)

## Synthetic demo design

The synthetic universe (`src/synthetic.rs`) generates four series of length N.
Two of them (`A_lead`, `A_lag`) share a Pareto-distributed upper-tail shock so
that χ for that pair tends to a non-trivial positive limit; the other two
(`B_indep1`, `B_indep2`) receive independent shocks. The bulk of every series
is Gaussian. A fixed seed (`SEED = 2026_04_25`) makes the run reproducible.

## Tests

```bash
cargo test --release
```

Runs an integration test that drives `evt-cli --synthetic` and asserts the
output JSON has the expected schema and finite GPD parameters per asset.

## References

- Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values.* Springer.
- Embrechts, P., Klüppelberg, C. & Mikosch, T. (1997). *Modelling Extremal Events.*
- Pickands, J. (1975). *Statistical inference using extreme order statistics.*
- Hosking, J. R. M. & Wallis, J. R. (1987). *Parameter and quantile estimation for the generalized Pareto distribution.*

## License

MIT © Daniel Vieira Gatto.
