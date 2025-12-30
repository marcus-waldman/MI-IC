# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **statistical methodology research project** focused on corrected information criteria for model selection under multiple imputation (MI). The project develops theoretical foundations and practical methods for properly penalizing model complexity when selecting among candidate models fitted to multiply-imputed datasets.

### Key Research Problem
Existing model selection methods (AIC, BIC) when applied to imputed data incur an **imputation bias** that scales with model complexity, favoring overfitting. The project derives corrected penalties that properly account for this bias.

### Directories
- **`claude/`**: Research documentation and paper drafts
  - `paper_outline_v4.md`: Formal research paper structure with full methodology, theory, and key results
  - `mi_deviance_bias_derivation.html`: Interactive technical document with step-by-step mathematical derivations
- **`simulations/`**: Directory for simulation code and results (currently empty; will contain R/Python simulation studies)

## Development Setup

### R Environment
- R version: 4.5.1 (from `C:\Program Files\R\R-4.5.1\bin`)
- Primary language for statistical computation and simulations
- Use `Rscript.exe` for batch execution

### Running Simulations
When simulation code is developed:
```bash
Rscript path/to/simulation.R
```

For interactive development:
```bash
R
```

## Project Architecture

### Theoretical Foundation
The core methodology is based on decomposing total bias when using multiply-imputed data for model selection:

```
Total Bias = Imputation Bias + Optimism
           = tr(RIV) + Q + tr(RIV)    [for AIC]
           = tr(RIV) + (Q log N)/2    [for BIC, approximation]
```

Where:
- **Imputation Bias**: Arises from double dipping (using data to estimate parameters, then generate imputations from those estimates)
- **Optimism**: Standard Efron-type bias from training/test split
- **RIV** = Relative Increase in Variance = `(1 + 1/M) * W^{-1} * B`
  - W = within-imputation variance
  - B = between-imputation variance
  - M = number of imputations

### Corrected Criteria
The key results are:
- **MI-AIC**: Penalty = `2Q + 4*tr(RIV)`
- **MI-BIC**: Penalty = `Q*log(N) + 2*tr(RIV)`

Both are **directly computable from standard MI output** (W and B matrices).

### Comparison to Existing Methods
The project shows that AICcd, PDIO, and Claeskels & Consentino's methods underpenalize by omitting the imputation bias correction term. See paper section 7 for detailed comparison.

## Mathematical Notation & Key Relationships

- `Q` = number of parameters in candidate model
- `N` = sample size
- `θ` = parameter vector
- `ℓ_com(θ)` = complete-data log-likelihood
- `ℓ_obs(θ)` = observed-data log-likelihood
- `Q(θ|θ̃)` = Q-function (expected complete-data log-likelihood)
- `ℐ_com`, `ℐ_obs`, `ℐ_mis|obs` = Fisher information matrices
- `W^{-1}` should be regularized if model is high-dimensional

## Integration with Paper Outline

The `paper_outline_v4.md` file contains:
1. **Abstract** summarizing the core problem and solution
2. **Sections 1-2**: Setup, notation, and foundational definitions
3. **Section 3**: Derivation of imputation bias formula
4. **Sections 4-5**: MI-AIC and MI-BIC derivations
5. **Section 6**: Extension to missing covariates (conditional criteria)
6. **Section 7**: Comparison to existing methods
7. **Section 8**: Simulation study design
8. **Section 9**: Discussion and practical guidance

The HTML document provides step-by-step derivations not in the outline.

## Future Work Hints

Based on Section 8, when developing simulation code, expect to:
- Generate synthetic data with nested model structures
- Implement MAR (Missing At Random) missingness mechanisms
- Compare selection replication rates across methods
- Test varying sample sizes (N ∈ {100, 500, 1000}) and missing rates (10%, 30%, 50%)
- Verify that MI-AIC/MI-BIC achieve higher replication of complete-data decisions than existing methods
