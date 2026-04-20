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

Requires R >= 4.0 with packages: `lavaan`, `mice`, `MASS`, `SeedMaker`,
`pbapply`, `parallel`.

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
| `compute_all_ic()` | Compute all 7 IC methods for one model |
| `load_results()` | Load saved results into a flat data.frame |
| `selection_accuracy_table()` | % selecting true model (M1) by condition |
| `plot_accuracy_by_n()` | Quick diagnostic plot |

## IC methods computed

1. **AIC_com / BIC_com** — Oracle (complete-data)
2. **AIC_adhoc / BIC_adhoc** — Mean of per-imputation IC
3. **AICcd** — Cavanaugh & Shumway (penalty: 2p + 2·tr(RIV))
4. **MI-AIC** — Proposed (penalty: 2p + tr(RIV), eval at pooled θ̄)
5. **MI-BIC** — Proposed (penalty: p·log(N) + tr(RIV), eval at pooled θ̄)

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
