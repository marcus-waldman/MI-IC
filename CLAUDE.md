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

### Critical Assumption: Congeniality (verified 2026-04-23)
The `+½tr(RIV)` bias result requires **congenial imputation**: the imputation model must contain the analysis model's parameter space. Empirically verified at N=250, mr=0.40, M=50:
- **Amelia (congenial, joint MVN)**: r(Term 1) = +0.85 ≈ theory's +1 (n=1989)
- **mice PMM (uncongenial, chained regressions)**: r(Term 1) = -0.39 (n=496, sign-flipped)

See `claude/notes/2026-04-23-congeniality-finding.md` for the full comparison, formal tests, and implications. The v4 derivation is correct under its stated congeniality assumption; uncongenial MI (standard mice PMM in SEM contexts) needs a separate correction.

### Manuscript Reframe: Chi-Square Recalibration (2026-04-26)
The bias correction's value-add is **NOT** model selection (where MI-AIC and AIC_adhoc give nearly identical % selecting M1, ~75% oracle agreement under both amelia and PMM). It is the **SEM goodness-of-fit chi-square test**.

At N=250, mr=0.40, M=50, true model M1 (df=22), Type I error at α=0.05:
- **chi2_adhoc**: 21% (4× nominal) — **standard practice**
- **chi2_MI**: 9.8% (with v4 first-moment correction)
- **chi2_MI,corr** (v4.5): pending validation, predicted ~5% (nominal)
- **chi2_D3** (Meng-Rubin): 14.8%
- **chi2_com** (oracle): 6.4%

MI bias correction fixes the **first moment** (E[chi2_MI] = 22.9 ≈ df=22, vs E[chi2_adhoc] = 27.9). Residual Type I inflation (10% vs 5%) is from **second-moment** inflation (Var[chi2_MI] = 63 vs Var[chi2_com] = 49) — closed in v4.5 §13 (see below).

See `claude/notes/2026-04-26-chisquare-recalibration.md`.

### v4.5 Closure: Finite-M Variance Correction (2026-04-29)
The v4.5 derivation closes the chi-square variance gap. The per-imputation chi-square decomposes as $\chi^2_m = A_m + B_m$ along the obs/mis split. As $M\to\infty$, $\bar A_\infty$ is the FIML observed-data chi-square (variance ≈ Var[χ²_com]) and $\bar B_\infty = 2\,\text{KL}(p_\text{imp} \,\|\, p_{M_1})$ — the imputation-flexibility KL divergence. Step 12 quadratic Taylor expansion + generalized symmetric eigendecomposition gives canonical Karhunen-Loève form $\bar B_\infty = \sum_{j=1}^{\Delta p}\lambda_j Z_j^2$ (generalized chi-square), where $\{\lambda_j\}$ are eigenvalues of $\text{RIV}_\perp = \mathcal{I}^\text{imp,mis}_{\perp\perp}(\mathcal{I}^\text{imp,obs}_{\perp\perp})^{-1}$.

**Closed form**: $\text{Var}[\chi^2_\text{MI}] = \text{Var}[\chi^2_\text{com}] + 2\sum_j\lambda_j^2 + O(M^{-1})$.

**Asparouhov-Muthén-style scaled-shifted correction**:
$$\chi^2_\text{MI,corr} = a\chi^2_\text{MI} + b, \quad a = \sqrt{\frac{2\,\text{df}}{2\,\text{df} + 2\sum_j\lambda_j^2}}, \quad b = \text{df}(1-a)$$

Implementation note: when computing the eigendecomposition of $W^{-1}B$ in R, symmetrize via Cholesky (`A_sym = solve(t(chol(W))) %*% B %*% solve(chol(W))`) and use `eigen(A_sym, symmetric = TRUE, only.values = TRUE)` for stability. See v4.5 §13 of `claude/derivations/mi_deviance_bias_derivation_v4.qmd`.

Open work: `miicsem` runtime implementation (`compute_chi2_MI_corrected()`) and end-to-end Type I validation. Plan: `~/.claude/plans/cozy-imagining-waterfall.md`.

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
- **Location**: `packages/miicsem/` (installable R package)
- **Design**: 3-factor CFA with 9 indicators
- **Sample sizes**: N = 100, 250, 500, 1000
- **Missing rates**: 10%, 25%, 40% MCAR
- **Imputations**: M = 100
- **Replications**: 2000 per condition
- **Methods compared**: Complete-data IC, Ad hoc, AICcd, D_LR-based AIC (Consentino & Claeskens 2010), MI-AIC, MI-BIC
- **Output**: Selection accuracy tables, IC distributions, figures; per-model deviances (DEV_com, DEV_adhoc, MI_DEVIANCE, MR_DEVIANCE) and three chi-squares (chi2_com, chi2_MI, chi2_D3) vs the saturated reference

### Study 3: LPA (Mixture Modeling)
- **Location**: `simulations/study3-lpa/`
- **Design**: 3-class LPA with 5-6 continuous indicators
- **Sample sizes**: N = 200, 500, 1000
- **Missing rates**: 15%, 30% MAR
- **Imputations**: M = 20
- **Replications**: 500 per condition
- **Methods compared**: Complete-data IC, Ad hoc, AICcd, D_LR-based AIC (Consentino & Claeskens 2010), MI-AIC, MI-BIC
- **Output**: Class enumeration accuracy, figures
