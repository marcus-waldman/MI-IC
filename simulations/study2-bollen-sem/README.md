# Study 2: SEM Simulation (Bollen et al. 2014)

Replicates Bollen et al. (2014) SIM1 design with missing data + multiple imputation to evaluate MI-AIC and MI-BIC against existing information criteria.

## Design

- **Model**: 3-factor SEM with 9 indicators and cross-loadings (Paxton et al., 2001)
- **12 competing models**: M1 (true), overspecified (M6,M8,M9), underspecified (M2-M4,M10), mixed (M5,M7), wrong structure (M11,M12)
- **Sample sizes**: N = 100, 250, 500, 1000, 5000
- **Missing rates**: 10%, 25%, 40% MCAR
- **Imputations**: M = 20 (PMM via mice)
- **Replications**: 1000 per condition (15 conditions total)

## IC Methods Compared

1. **AIC_com / BIC_com** — Oracle (complete-data)
2. **AIC_adhoc / BIC_adhoc** — Mean of per-imputation IC
3. **AICcd** — Cavanaugh & Shumway (penalty: 2p + 2·tr(RIV))
4. **MI-AIC** — Our proposal (penalty: 2p + tr(RIV), evaluated at pooled θ̄)
5. **MI-BIC** — Our proposal (penalty: p·log(N) + tr(RIV), evaluated at pooled θ̄)

## Usage

```r
# Quick validation (1 rep)
source("tests/test_single_rep.R")

# Full simulation
source("R/06_run_simulation.R")

# Analyze results
source("R/07_analyze_results.R")
```

## File Structure

- `R/00_config.R` — Design parameters, population covariance
- `R/01_model_specifications.R` — 12 lavaan model syntaxes
- `R/02_generate_and_ampute.R` — Data generation and MCAR amputation
- `R/03_fit_and_pool.R` — Core fitting, Rubin's rules pooling, fixed-param loglik
- `R/04_compute_ic.R` — 5 IC methods
- `R/05_run_replication.R` — Single-rep wrapper
- `R/06_run_simulation.R` — Parallel driver (pbapply + parallel)
- `R/07_analyze_results.R` — Tables and figures
- `tests/test_single_rep.R` — Quick validation test
