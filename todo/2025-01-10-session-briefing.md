# Session Briefing: January 10, 2025

## Where We Left Off

### 1. Three-Term Bias Decomposition (NEW)

Added Section 7 to `claude/derivations/mi_deviance_bias_derivation_v3_3.qmd` deriving the **practical bias** when practitioners evaluate log-likelihoods at imputation-specific MLEs θ̂^(m) rather than at the pooled estimate θ̂_obs.

**Key result:**

| Implementation | Bias | Penalty |
|----------------|------|---------|
| Theoretical (evaluate at θ̂_obs) | ½tr(RIV) | 2Q + tr(RIV) |
| Practical (evaluate at θ̂^(m)) | tr(RIV) | 2Q + 2tr(RIV) |

The three-term decomposition:
- **Term 0** = ½tr(RIV) — evaluation point variability
- **Term 1** = tr(RIV) — imputation bias
- **Term 2** = -½tr(RIV) — estimation mismatch

This explains why AICcd has penalty 2Q + 2tr(RIV) — it targets the practical implementation.

### 2. Proper vs Improper MI Simulation Code

Updated `simulations/study1-pedagogy-clean/` to compute both:
- **Improper MI**: Q(θ̂_obs | θ̃ = θ̂_obs) — θ̃ fixed at MLE
- **Proper Bayesian MI**: E_θ̃[Q(θ̂_obs | θ̃)] via delta method (Hessian correction)

Files modified:
- `R/02_ampute_data.R` — Simplified to monotone-only (MCAR, MAR, MNAR)
- `R/03_analytical_functions.R` — Added `compute_Q_analytical_proper()` with Hessian corrections
- `R/05_run_simulation.R` — Now returns both improper and proper results
- `test_clean.R` — Compares improper vs proper validation

---

## What Still Needs To Be Done

### Run the Simulation

```bash
cd simulations/study1-pedagogy-clean
Rscript test_clean.R
```

This will run 5000 replications validating that:
- Term 1 (improper) ≈ tr(RIV)
- Term 2 ≈ -½tr(RIV)
- Total (improper) ≈ ½tr(RIV)
- Term 1 (proper) ≈ tr(RIV)
- Total (proper) ≈ ½tr(RIV)
- Hessian correction is negative

### Important Note

**The current simulation validates the THEORETICAL bias (½tr(RIV)), not the PRACTICAL bias (tr(RIV)).**

The simulation evaluates Q at a single point θ̂_obs (the FIML estimate), which corresponds to our MI-AIC derivation. It does NOT simulate the practical workflow where each imputation's log-likelihood is evaluated at its own MLE θ̂^(m).

To validate the practical bias (Term 0 + Term 1 + Term 2 = tr(RIV)), you would need a separate simulation that:
1. Generates M imputed datasets
2. Fits the model to each to get θ̂^(m)
3. Computes ℓ^(m)(θ̂^(m)) for each
4. Averages and compares to ℓ_com(θ̂_com)

This is a potential future extension but not currently implemented.

---

## File Locations

- Derivation document: `claude/derivations/mi_deviance_bias_derivation_v3_3.qmd`
- Simulation code: `simulations/study1-pedagogy-clean/`
- Plan file: `C:\Users\marcu\.claude\plans\delightful-crafting-bear.md`
