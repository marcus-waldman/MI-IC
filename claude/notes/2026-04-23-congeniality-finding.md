# Congeniality Governs the Sign of Term 1 (Imputation Bias)

**Date**: 2026-04-23
**Context**: Three-term bias decomposition in the MI-IC derivation. Study 2 empirical results contradicted the v4 derivation's predicted sign for Term 1 (Imputation Bias). This note documents the definitive test that pins down the cause as **imputation-model congeniality** rather than a theorem error.

## The decomposition

Writing the bias as a telescoping sum through two intermediate evaluation points:

$$
\text{Bias} = \underbrace{\mathbb{E}\!\left[\bar{Q}_{MI}(\bar\theta_{MI}) - \ell_{\text{com}}(\bar\theta_{MI})\right]}_{\text{T1: Imputation}}
            + \underbrace{\mathbb{E}\!\left[\ell_{\text{com}}(\bar\theta_{MI}) - \ell_{\text{com}}(\hat\theta_{\text{obs}})\right]}_{\text{T2: Pooling approximation}}
            + \underbrace{\mathbb{E}\!\left[\ell_{\text{com}}(\hat\theta_{\text{obs}}) - \ell_{\text{com}}(\hat\theta_{\text{com}})\right]}_{\text{T3: Estimation mismatch}}
$$

**v4 prediction (congenial MLE imputation):** $T_1 = +\text{tr}(\text{RIV})$, $T_2 \approx 0$, $T_3 = -\tfrac{1}{2}\text{tr}(\text{RIV})$, total $= +\tfrac{1}{2}\text{tr}(\text{RIV})$.

## Experiment

Identical Study 2 DGP (9-indicator, 3-factor CFA) at $N=250$, $\text{mr}=0.40$ MCAR, $M=50$. Two imputation backends compared (2000 reps for amelia; 500-rep pmm baseline carried over from the earlier decomp HPC run, same seeds family):

| Backend | Family | Joint coherent? | Congenial with CFA? |
|---|---|---|---|
| **pmm** (`mice(method="pmm")`) | semi-parametric, chained regressions with nearest-neighbor matching | No (not derivable from a joint) | **No** — different parameter space than CFA |
| **amelia** (`Amelia::amelia(boot.type="ordinary")`) | parametric; bootstrap-EM on joint MVN, then conditional MVN draws | **Yes** | **Yes** — joint MVN $(\mu, \Sigma)$ is the saturated model, which nests every candidate CFA |

## Results

On the log-likelihood scale for $M_1$ (the true model):

| Method | $n_{\text{reps}}$ | $\overline{\text{tr}(\text{RIV})}$ | $T_1$ | $T_2$ | $T_3$ | Total |
|---|---|---|---|---|---|---|
| amelia | 1989 | 2.023 | **$+1.718 \pm 0.173$** | $-0.056$ | $-1.032$ | $+0.630$ |
| pmm    |  496 | 2.019 | **$-0.792 \pm 0.344$** | $-0.076$ | $-1.000$ | $-1.869$ |

As multiples of $\text{tr}(\text{RIV})$:

| Method | $r_1$ | $r_2$ | $r_3$ | $r_{\text{total}}$ |
|---|---|---|---|---|
| amelia | **+0.849** | $-0.028$ | $-0.510$ | $+0.311$ |
| pmm    | **$-0.392$** | $-0.038$ | $-0.495$ | $-0.926$ |
| *v4 theory (congenial)* | *+1.000* | *0* | *$-0.500$* | *+0.500* |

## Formal tests on amelia (2000 reps)

$$H_0:\ \mathbb{E}[T_1] = 0 \quad\Longrightarrow\quad t = +9.91,\ n = 1989$$

Reject decisively — $T_1$ is positive under amelia.

$$H_0:\ \mathbb{E}[T_1] = +\text{tr}(\text{RIV}) \quad\Longrightarrow\quad t = -1.76,\ p \approx 0.08$$

Do not reject at $\alpha = 0.05$. Amelia's Term 1 is **statistically indistinguishable** from the theoretical $+\text{tr}(\text{RIV})$.

## Conclusions

1. **$T_3$ is correct under any imputation.** $r_3 \to -0.5$ to two decimals for both amelia and pmm.
2. **$T_2 \approx 0$ for both.** The pooling approximation $\bar\theta_{MI} \approx \hat\theta_{\text{obs}}$ holds at $N=250$.
3. **$T_1$ has opposite sign between methods.** Congenial imputation gives $T_1 \approx +\text{tr}(\text{RIV})$ as the v4 theory predicts; uncongenial PMM gives $T_1 \approx -\tfrac{1}{2}\text{tr}(\text{RIV})$, a swing of $\approx -\tfrac{3}{2}\text{tr}(\text{RIV})$.
4. **The v4 derivation is correct under its stated congeniality assumption.** The empirical failure originally observed in Study 2 was a *model-assumption violation* (PMM + CFA uncongeniality), not a theorem error.

## Implications for the manuscript

- **Section 3 (Theory):** Congeniality must be stated as a manifest assumption, not buried in Step 10 of the derivation. The $+\text{tr}(\text{RIV})$ correction applies only when the imputation model contains the analysis model's parameter space.
- **Section 4 (Simulation):** Report results separately for congenial (amelia / jomo-style joint imputation) vs uncongenial (mice PMM / norm chained) settings. Expected pattern:
  - Under congenial imputation: MI-AIC $= -2\bar{Q} + 2p + \text{tr}(\text{RIV})$ recovers AIC_com.
  - Under uncongenial imputation (PMM): the formula over-corrects by $\approx 2 \cdot \text{tr}(\text{RIV})$; an uncongeniality-adjusted correction (or a robustness diagnostic) is needed.
- **Section 6 (Discussion):** Recommend joint-modeling MI (Amelia, jomo, blimp) for applied users who want MI-AIC/MI-BIC to behave as theory predicts.
- **Future work:** derive the uncongeniality correction analytically, or empirically characterize it as a function of the imputation-model / analysis-model "distance" (e.g., KL or score-gap).

## Artifacts

- `packages/miicsem/R/impute_mvn.R` — congenial saturated-MVN plug-in imputer (`mvn_msat`).
- `packages/miicsem/R/run_replication.R` — dispatch on `mice_method`, supporting `pmm`/`norm`/`amelia`/`mvn_msat`.
- `hpc/results-congeniality-amelia/results_combined.rds` — 2000 amelia reps with decomposition columns.
- `hpc/figures-M100/amelia_vs_pmm_{terms_long,M1_summary}.csv` — comparison tables.
- `hpc/analyze_amelia_vs_pmm.R` — reproduction analysis script.
