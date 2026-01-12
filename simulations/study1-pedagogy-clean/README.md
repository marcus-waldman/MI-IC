# Study 1: Pedagogical Validation of MI-AIC Bias Formula

**Goal**: Validate that the bias formula **Bias ≈ ½tr(RIV)** holds for multiply-imputed data under different imputation model specifications.

## Overview

This simulation validates the theoretical derivation that the Q-function (pooled MI deviance) has a bias of **½tr(RIV)** when evaluated at the observed-data MLE, where RIV is the Relative Increase in Variance from multiple imputation.

**Key Features:**
- **Analytical approach** (M = ∞): Exact computation via Missing Information Principle, no MCMC sampling
- **Three imputation models**: Tests both congenial (Models A, C) and uncongenial (Model B) scenarios
- **Proper Bayesian MI**: Includes Hessian corrections for integration over posterior distribution
- **Fast and exact**: Validated against brms Monte Carlo integration (z-scores << 2.0)

---

## Quick Start

### Option 1: Quick Test (1 replication, ~5 seconds)

Tests all three models with a single dataset to verify everything works:

```r
source("test_single_rep.R")
```

### Option 2: Small Simulation (100 reps, ~1-2 minutes)

```r
source("R/00_seed_management.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/03_analytical_functions.R")
source("R/04_compute_metrics.R")
source("R/05_run_simulation.R")

seeds <- generate_seeds(master_seed = 12345, n_reps = 100)

results_A <- run_simulation(
  seeds = seeds,
  n = 500,
  mechanism = "MCAR",
  gamma = 0.5,
  cut1 = 0.5,
  cut2 = 0.75,
  true_params = list(a = 0.5, b = 0.5, sigma2_M = 1, sigma2_Y = 1),
  imputation_model = "A",
  n_cores = 4
)

# Check validation
mean(results_A$term1_improper) / mean(results_A$tr_RIV)  # Should be ~1.0
```

### Option 3: Full Validation (5000 reps, ~15-20 minutes)

```r
source("test_clean.R")
```

---

## Three Imputation Models

| Model | Structure | Parameters | Congeniality | Expected | Use Case |
|-------|-----------|------------|--------------|----------|----------|
| **A** | M\|X, Y\|M (matches DGP) | (a, b, σ²_M, σ²_Y) | ✓ Congenial | Bias ≈ ½tr(RIV) | Validate formula under correct specification |
| **B** | M\|X, Y marginal (ignores M→Y) | (a, σ²_M, μ_Y, τ²_Y) | ✗ Uncongenial | Bias ≠ ½tr(RIV) | Demonstrate uncongeniality breaks formula |
| **C** | (X,M,Y) ~ N₃(μ, Σ) saturated | 9 params (3 means + 6 covariances) | ✓ Asymptotically congenial | Bias ≈ ½tr(RIV) | Validate formula with flexible specification |

**Key Insight:** Model B's uncongeniality manifests as ignoring the observed mediator M when imputing Y, and incorrectly assuming Cov[M,Y|X] = 0.

---

## Installation

### Requirements

- **R version**: 4.0.0 or higher (tested on 4.5.1)
- **Required packages**:
  ```r
  install.packages(c("pbapply", "parallel"))
  ```
- **Optional** (for validation against brms):
  ```r
  install.packages(c("brms", "posterior"))
  ```

---

## Usage

### Running Test Scripts

**Quick test (1 rep):**
```r
source("test_single_rep.R")
```
Output: Comparison table showing all three models, expected runtime ~5 seconds

**Full validation (5000 reps):**
```r
source("test_clean.R")
```
Output: Validation results with pass/fail for each model, expected runtime ~15-20 minutes

**brms validation:**
```r
source("test_analytical_integration.R")
```
Output: Comparison of analytical vs Monte Carlo integration, expected runtime ~10 minutes

### Custom Simulation

```r
# Source all functions
source("R/00_seed_management.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/03_analytical_functions.R")
source("R/04_compute_metrics.R")
source("R/05_run_simulation.R")

# Configure
n_reps <- 1000
n <- 500
mechanism <- "MCAR"  # or "MAR" or "MNAR"
model <- "A"         # or "B" or "C"

# Generate reproducible seeds
seeds <- generate_seeds(master_seed = 12345, n_reps = n_reps)

# Run simulation
results <- run_simulation(
  seeds = seeds,
  n = n,
  mechanism = mechanism,
  gamma = 0.5,       # Strength parameter for MAR/MNAR
  cut1 = 0.5,        # Quantile for Y-only missing
  cut2 = 0.75,       # Quantile for M+Y missing
  true_params = list(a = 0.5, b = 0.5, sigma2_M = 1, sigma2_Y = 1),
  imputation_model = model,
  n_cores = parallel::detectCores() - 1
)

# Validate
cat(sprintf("Mean term1/tr(RIV): %.3f (should be ~1.0 for A,C)\n",
            mean(results$term1_improper) / mean(results$tr_RIV)))
```

### Expected Runtime

| Replications | Model A | Model B | Model C | Total (all 3) |
|--------------|---------|---------|---------|---------------|
| 1 | <1 sec | <1 sec | <1 sec | ~5 sec |
| 100 | ~20 sec | ~20 sec | ~20 sec | ~1 min |
| 1000 | ~3 min | ~3 min | ~3 min | ~10 min |
| 5000 | ~15 min | ~15 min | ~15 min | ~45 min |

*Tested on 8-core machine with parallel execution*

---

## Configuration

### Sample Size
```r
n <- 500  # Change to 100, 250, 1000, etc.
```

### Missingness Mechanism
```r
mechanism <- "MCAR"  # Missing Completely At Random
mechanism <- "MAR"   # Missing At Random (depends on X)
mechanism <- "MNAR"  # Missing Not At Random (depends on M, Y)
```

### Missingness Rates
```r
cut1 <- 0.5   # Quantile for Y-only missing (higher = more missing)
cut2 <- 0.75  # Quantile for M+Y missing (higher = more missing)
gamma <- 0.5  # Strength of MAR/MNAR relationship (0 = MCAR, 1 = strong)
```

### True Parameters
```r
true_params <- list(
  a = 0.5,         # M <- X slope
  b = 0.5,         # Y <- M slope
  sigma2_M = 1,    # Residual variance for M
  sigma2_Y = 1     # Residual variance for Y
)
```

---

## Output Interpretation

### Key Columns

| Column | Meaning |
|--------|---------|
| `term1_improper` | Imputation bias: E[Q̄_MI(θ̂_obs)] - E[ℓ_com(θ̂_obs)] |
| `term2` | Estimation mismatch: E[ℓ_com(θ̂_obs)] - E[ℓ_com(θ̂_com)] |
| `total_improper` | Total bias (improper MI): term1 + term2 |
| `hessian_correction` | Correction from integrating over posterior |
| `total_proper` | Total bias (proper MI): includes Hessian correction |
| `tr_RIV` | Trace of Relative Increase in Variance matrix |
| `Q_improper` | Q-function evaluated at point estimate |
| `Q_proper` | Q-function with posterior integration |

### Expected Patterns

**Model A (Congenial):**
```
mean(term1_improper) / mean(tr_RIV) ≈ 1.0
mean(total_improper) / (0.5 * mean(tr_RIV)) ≈ 1.0
```

**Model B (Uncongenial):**
```
mean(term1_improper) / mean(tr_RIV) ≠ 1.0  (demonstrates uncongeniality)
```

**Model C (Saturated MVN):**
```
mean(term1_improper) / mean(tr_RIV) ≈ 1.0
mean(total_improper) / (0.5 * mean(tr_RIV)) ≈ 1.0
```

**All Models:**
```
mean(term2) / (-0.5 * mean(tr_RIV)) ≈ 1.0
mean(hessian_correction) < 0  (log-likelihood concavity)
```

---

## Directory Structure

```
study1-pedagogy-clean/
├── R/
│   ├── 00_seed_management.R          # SeedMaker for reproducible seeds
│   ├── 01_generate_data.R             # Data generation from mediation model
│   ├── 02_ampute_data.R               # Create MCAR/MAR/MNAR missingness
│   ├── 03_analytical_functions.R      # Core: Q-function, Hessian corrections
│   ├── 04_compute_metrics.R           # Metrics computation
│   └── 05_run_simulation.R            # Main loop with parallel execution
├── test_single_rep.R                  # Quick test (1 rep, all models)
├── test_clean.R                       # Full validation (5000 reps)
├── test_analytical_integration.R      # brms MC validation
├── IMPLEMENTATION_SUMMARY.md          # Technical documentation
└── README.md                          # This file
```

---

## Troubleshooting

**Error: "cannot open the connection"**
- Check that you're in the correct directory: `setwd("simulations/study1-pedagogy-clean")`

**Simulation runs slowly**
- Reduce `n_cores` if memory constrained
- Use smaller `n` or fewer reps for testing
- Close other applications

**Ratios not close to 1.0**
- With small replications (< 100), expect sampling variability
- Use at least 1000 reps for stable validation
- Check that model specification matches expectation

**brms validation fails**
- Install brms: `install.packages("brms")`
- Requires Stan compiler (can take 10+ minutes on first run)

---

## Technical Details

For implementation details, mathematical derivations, and validation results, see:
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Complete technical documentation
- **Manuscript**: `../../manuscript/mi-ic-article.qmd`
- **Derivation**: `../../claude/derivations/mi_deviance_bias_derivation_v3_3.qmd`

---

## References

**Key Result:**
```
Total Bias = E[Q̄_MI(θ̂_obs)] - E[ℓ_com(θ̂_com)]
           = tr(RIV) + (-½tr(RIV))
           = ½tr(RIV)
```

where:
- **Term 1**: tr(RIV) — Imputation bias
- **Term 2**: -½tr(RIV) — Estimation mismatch
- **RIV** = Relative Increase in Variance = (1 + 1/M) W⁻¹B

**Corrected Information Criteria:**
- MI-AIC: Penalty = 2Q + tr(RIV)
- MI-BIC: Penalty = Q·log(N) + tr(RIV)
