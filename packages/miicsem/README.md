# miicsem

Simulation framework for evaluating **MI-AIC** and **MI-BIC** information
criteria against **AICcd** and ad hoc alternatives in structural equation
models with missing data and multiple imputation.

Implements the Bollen et al. (2014) SIM1 design (3-factor SEM with
cross-loadings, 9 indicators, 12 competing models). Designed for
reproducible execution on high-performance computing (HPC) environments.

## Installation

```r
# From GitHub (development version)
devtools::install_github("marcus-waldman/MI-IC", subdir = "packages/miicsem")
```

Requires R >= 4.0. Runtime dependencies (`lavaan`, `mice`, `MASS`,
`SeedMaker`, `pbapply`) are installed automatically on first
`run_simulation()` call if missing — opt out with
`run_simulation(install_missing = FALSE)`.

## Quick start

```r
library(miicsem)

# Pilot run (10 reps, single condition, ~3 min on 16 cores)
pilot <- run_simulation(
  n_reps       = 10,
  sample_sizes = 250,
  miss_rates   = 0.25
)

# Full simulation (1000 reps x 15 conditions)
results <- run_simulation(
  n_reps      = 1000,
  n_cores     = parallel::detectCores(),
  seed        = 32897891,
  results_dir = "sim_results/"
)
```

## Function reference

| Function | Purpose |
|----------|---------|
| `run_simulation()` | Main entry point; runs full grid with checkpointing |
| `run_one_rep()` | Single replication (useful for debugging) |
| `get_config()` | Retrieve default configuration |
| `get_sim1_models()` | 12 lavaan model syntaxes |
| `get_saturated_model()` | Saturated (unstructured) lavaan syntax |
| `compute_pop_starts()` | Pre-compute warm-start coefficients from Σ_pop |
| `compute_all_ic()` | Compute all 7 IC methods for one model |
| `compute_D3()` | Meng-Rubin (1992) D_3 for one nested model pair |
| `load_results()` | Selection winners per rep |
| `load_deviances()` | Long-form deviance table |
| `load_chi_squares()` | Long-form chi-square table |
| `chi2_summary_table()` | Mean chi-squares by (N, miss_rate, model) |
| `selection_accuracy_table()` | % selecting true model (M1) by condition |
| `plot_accuracy_by_n()` | Quick diagnostic plot |

## IC methods computed

1. **AIC_com / BIC_com** — Oracle (complete-data)
2. **AIC_adhoc / BIC_adhoc** — Mean of per-imputation IC
3. **AICcd** — Cavanaugh & Shumway (penalty: 2p + 2·tr(RIV))
4. **MI-AIC** — Proposed (penalty: 2p + tr(RIV), eval at pooled θ̄)
5. **MI-BIC** — Proposed (penalty: p·log(N) + tr(RIV), eval at pooled θ̄)

## Deviances saved per model

All four on the deviance (−2 ℓ) scale. Computed for all 12 candidates
and the saturated model:

| Quantity | Formula |
|----------|---------|
| `DEV_com` | −2 ℓ_com(θ̂_com) — oracle |
| `DEV_adhoc` | mean over imputations of −2 ℓ_m(θ̂_m) — naive |
| `MI_DEVIANCE` | −2 Q̄(θ̄) + tr(RIV) — our debiased |
| `MR_DEVIANCE` | Meng-Rubin, anchored at saturated |

## Chi-squares saved per candidate model

All computed against the saturated model:

| Quantity | Formula | Role |
|----------|---------|------|
| `chi2_com` | DEV_com(M_j) − DEV_com(Msat) | Oracle chi-square |
| `chi2_MI`  | MI_DEVIANCE(M_j) − MI_DEVIANCE(Msat) | Ours; multivariate tr(RIV) |
| `chi2_D3`  | d̃̄ / (1 + r_L) | Meng-Rubin (1992); scalar r_L |
| `df`       | npar(Msat) − npar(M_j) | Same for all three |

## Warm starts

Each simulation run precomputes "infinite-N projection" parameter values
by fitting every candidate + saturated model once to `sigma_pop` with
`sample.nobs = 1e6`, then passes those coefficients as starting values
to every per-imputation fit. Purely a performance optimization; the MLEs
are unchanged.

## HPC usage

```bash
#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load R/4.5.1

Rscript -e 'miicsem::run_simulation(n_reps = 1000, results_dir = "/scratch/$USER/sim/")'
```

## Reproducibility

All randomness is seeded via `SeedMaker` (three independent streams per
replication: data generation, amputation, imputation). Default master
seed is `32897891`.

## Citation

See the parent repository [MI-IC](https://github.com/marcus-waldman/MI-IC)
for the underlying theoretical derivation and manuscript.
