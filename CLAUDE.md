# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **statistical methodology research project** developing a principled derivation of information criteria (AIC, BIC) for multiply-imputed data. The goal is **complete-data replication**—deriving IC that answers "what would I have concluded had my data been complete?"

### Key Research Contribution
The project provides a theoretically sound derivation showing that the Q-function (pooled deviance across imputations) is biased for the complete-data log-likelihood at the complete-data MLE. This bias equals ½tr(RIV), leading to corrected penalties that better replicate complete-data inference.

### Directory Structure (Post-Reorganization)

```
MI-IC/
├── manuscript/                          # Article in production
│   ├── mi-ic-article.qmd               # Main manuscript file
│   ├── references.bib                  # Bibliography (linked to Zotero)
│   ├── figures/                        # Publication-ready figures
│   └── tables/                         # Publication-ready tables
├── simulations/                        # Simulation code
│   ├── study2-bollen-sem/              # Study 2: SEM with missing data + MI
│   ├── study3-lpa/                     # Study 3: Latent Profile Analysis
│   ├── premise-tests/                  # Archived premise/verification simulations
│   │   ├── 000-premise-very-basic-Y/
│   │   └── 00-premise-basic-Y/
│   └── utils/                          # Shared R utilities
├── empirical-example/                  # Real data analysis
├── supplementary/                      # Supplementary materials (if needed)
├── claude/                             # Working drafts & development
│   ├── derivations/                    # Mathematical derivations
│   │   ├── mi_deviance_bias_derivation_v4.qmd    # Current version (v4.1)
│   │   ├── mi_deviance_bias_derivation_v4.html
│   │   ├── mi_deviance_bias_derivation_v4_files/
│   │   └── archive/                    # Older versions (v2, v2_1)
│   ├── outlines/                       # Article outlines
│   │   └── article_outline.qmd         # Detailed article structure plan
│   └── notes/                          # Other working notes
├── submission/                         # Submission materials
└── CLAUDE.md                           # This file
```

**Key File Locations:**
- **Main manuscript**: `manuscript/mi-ic-article.qmd`
- **Current derivation**: `claude/derivations/mi_deviance_bias_derivation_v4.qmd`
- **Article outline**: `claude/outlines/article_outline.qmd`
- **Bibliography**: `manuscript/references.bib` (link to Zotero MI-IC collection)

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
The target is **AIC_com**—the AIC that would have been computed at the complete-data MLE θ̂_com. The bias of the Q-function (pooled MI deviance) decomposes as:

```
Total Bias = E[Q̄_MI(θ̂_obs)] - E[ℓ_com(θ̂_com)]
           = Term 1 (Imputation Bias) + Term 2 (Estimation Mismatch)
           = tr(RIV) + (-½tr(RIV))
           = ½tr(RIV)
```

Where:
- **Term 1 (Imputation Bias)**: tr(RIV) — arises from imputing with estimated parameters
- **Term 2 (Estimation Mismatch)**: -½tr(RIV) — loss from evaluating at θ̂_obs instead of θ̂_com
- **RIV** = Relative Increase in Variance = `(1 + 1/M) * W^{-1} * B`

### Corrected Criteria
The bias correction applies to the deviance, not the penalty structure:
- **MI-AIC**: Penalty = `2Q + tr(RIV)`
- **MI-BIC**: Penalty = `Q*log(N) + tr(RIV)`

Both are **directly computable from standard MI output** (W and B matrices).

### Comparison to Existing Methods
AICcd (Cavanaugh & Shumway) has penalty `2Q + 2tr(RIV)`. The difference arises because AICcd targets predictive performance at θ̂_obs, not replication of complete-data inference at θ̂_com.

### Scope Decision: D_LR Connection
Consentino & Claeskens (2010) built an AIC from Meng & Rubin's (1992) D_L statistic: `aic(S, S_0) = -D_S + 2p_S`. This is a direct comparator that uses D_L's scalar r_L correction vs. our multivariate tr(RIV). **Decision**: Include D_LR-based AIC as a simulation comparator (Section 4.4). Introduce D_L in background (Section 2.4). Compare correction approaches in theory (Section 3.7). Discuss empirical results in Section 6.2. Flag corrected LR *tests* (not IC) as future work in Section 6.5.

## Mathematical Notation & Key Relationships

- `Q` = number of parameters in candidate model
- `N` = sample size
- `θ` = parameter vector
- `ℓ_com(θ)` = complete-data log-likelihood
- `ℓ_obs(θ)` = observed-data log-likelihood
- `Q(θ|θ̃)` = Q-function (expected complete-data log-likelihood)
- `ℐ_com`, `ℐ_obs`, `ℐ_mis|obs` = Fisher information matrices
- `W^{-1}` should be regularized if model is high-dimensional

## Target Journals
- *Structural Equation Modeling: A Multidisciplinary Journal*
- *Psychological Methods*

## Article Development Workflow

### Phase 0: File Organization ✓ (COMPLETE)
- Established clean directory structure separating article production from development
- Moved outline and derivations to working directories
- Created manuscript placeholder with YAML header and section stubs
- Archived older derivation versions

### Phase 1: Theory Development (IN PROGRESS)
1. **Verify derivation**: Ensure `claude/derivations/mi_deviance_bias_derivation_v4.qmd` is complete and rigorous
2. **Draft introduction**: Section 1 of `manuscript/mi-ic-article.qmd`
3. **Draft background**: Section 2 of manuscript
4. **Draft theory**: Section 3 of manuscript (use v3_3.qmd derivation as basis)

### Phase 2: Empirical Validation (PENDING)
5. **Study 2 (SEM)**: Implement in `simulations/study2-bollen-sem/`
   - Replicate Bollen et al. (2014) design with missing data + MI
   - Compute all IC variants (complete-data, ad hoc, MI-AIC, MI-BIC, AICcd)
   - Generate results tables and figures

6. **Study 3 (LPA)**: Implement in `simulations/study3-lpa/`
   - Design mixture model simulation with specified parameters
   - Compute all IC variants
   - Generate results tables and figures

7. **Empirical Example**: Implement in `empirical-example/`
   - Identify public behavioral science dataset
   - Apply all IC methods
   - Create results tables for Section 5 of manuscript

### Phase 3: Writing & Finalization (PENDING)
8. **Write simulations section**: Section 4 of manuscript (with embedded results)
9. **Write empirical example section**: Section 5 of manuscript
10. **Write discussion**: Section 6 of manuscript
11. **Finalize references**: Export from Zotero MI-IC collection to `manuscript/references.bib`
12. **Create publication figures**: Place in `manuscript/figures/`

## Simulation Studies

### Study 2: SEM (Replicating Bollen et al. 2014)
- **Location**: `simulations/study2-bollen-sem/`
- **Design**: 3-factor CFA with 9 indicators
- **Sample sizes**: N = 100, 250, 500, 1000, 5000
- **Missing rates**: 10%, 25%, 40% MCAR
- **Imputations**: M = 20
- **Replications**: 1000 per condition
- **Methods compared**: Complete-data IC, Ad hoc, AICcd, D_LR-based AIC (Consentino & Claeskens 2010), MI-AIC, MI-BIC
- **Output**: Selection accuracy tables, IC distributions, figures

### Study 3: LPA (Mixture Modeling)
- **Location**: `simulations/study3-lpa/`
- **Design**: 3-class LPA with 5-6 continuous indicators
- **Sample sizes**: N = 200, 500, 1000
- **Missing rates**: 15%, 30% MAR
- **Imputations**: M = 20
- **Replications**: 500 per condition
- **Methods compared**: Complete-data IC, Ad hoc, AICcd, D_LR-based AIC (Consentino & Claeskens 2010), MI-AIC, MI-BIC
- **Output**: Class enumeration accuracy, figures
