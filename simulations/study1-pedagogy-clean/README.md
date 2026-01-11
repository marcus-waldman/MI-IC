# Study 1: Pedagogical Validation of MI-AIC/MI-BIC Bias Formula

**Goal**: Validate that Bias ≈ ½tr(RIV) for multiply-imputed data

## Analytical Approach (M = Inf)

This implementation uses exact analytical computation of the Q-function and tr(RIV) via the Missing Information Principle, avoiding Monte Carlo sampling entirely.

### Key Features

- **Fast**: No MCMC sampling required
- **Exact**: Analytical computation of Q-function
- **Validated**: Matches theoretical predictions (ratios ≈ 1.0)

### Directory Structure

```
study1-pedagogy-clean/
├── R/
│   ├── 01_generate_data.R          # Data generation
│   ├── 02_ampute_data.R             # Create missingness
│   ├── 03_analytical_functions.R    # Core analytical functions
│   ├── 04_compute_metrics.R         # Metrics computation
│   └── 05_run_simulation.R          # Main simulation loop
├── config.R                         # Simulation configuration
├── run_study1.R                     # Main entry point
└── README.md                        # This file
```

## Usage

```r
# Run simulation
source("run_study1.R")
```

## Configuration

Set `M = Inf` in `config.R` to use analytical approach.
