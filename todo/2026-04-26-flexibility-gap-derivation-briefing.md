# Session Briefing: Derive the Flexibility-Gap Correction

**For**: a fresh-context Claude Code session
**Created**: 2026-04-26
**Goal**: derive the analytic form of the third finite-$N$ correction to the MI chi-square statistic — the **imputation-vs-analysis flexibility-gap correction** ($\Delta$ in v4.2's Section 11)

---

## TL;DR

Under standard MI workflows in SEM (Amelia/jomo's saturated joint MVN imputation, mice PMM's chained equations), the chi-square statistic $\chi^2_M$ for a candidate analysis model $M$ against the saturated reference $M_{\text{sat}}$ is **inflated** by an amount that

- is **not** explained by Bartlett-type ML finite-$N$ bias,
- is **not** explained by Wishart-type inverse-Hessian bias,
- empirically equals $\approx +5$ chi-square units at $N=250$, $q=9$ indicators, $p_{\text{analysis}}=32$, $p_{\text{sat}}=54$, mr=0.40,
- vanishes when the imputation model is constrained to the analysis model's parameter space (`mvn_M1` true-model imputation).

Your job: **derive the closed-form expression for this inflation as a function of $(p_{\text{imp}}, p_{\text{analysis}}, \text{mr}, N)$**.

---

## State at handoff

### What is already done

1. **v4 derivation** (`claude/derivations/mi_deviance_bias_derivation_v4.qmd`, version 4.2): single-model deviance bias $\mathbb{E}[\bar Q_{MI}] - \mathbb{E}[\ell_{\text{com}}] = +\tfrac{1}{2}\text{tr}(\text{RIV})$, plus a new Section 10 covering Bartlett and Wishart corrections, and a new Section 11 documenting the flexibility gap as future work (your work).

2. **Empirical verification** at $N=250$, mr=0.40, $M=50$, true model $M_1$ (df=22), $\alpha=0.05$:

   | Imputer | $\chi^2_{\text{com}}$ | $\chi^2_{\text{adhoc}}$ | $\chi^2_{\text{MI}}$ raw | $\chi^2_{\text{MI}}\cdot W$ | $\chi^2_{\text{MI}}\cdot W\cdot B$ |
   |---|---|---|---|---|---|
   | amelia (saturated) | 6.4% | 21.8% | 9.8% | 16.6% | 14.0% |
   | PMM (chained) | 7.1% | 20.8% | 11.1% | 16.3% | 14.5% |
   | **mvn_M1 (true model)** | 7.4% | **5.8%** | 2.2% | **4.8%** | 3.6% |

   Under truly congenial `mvn_M1` imputation, Wishart-corrected $\chi^2_{\text{MI}}$ achieves nominal Type I (4.8%). Under saturated/PMM imputation, NONE of {raw, Wishart, Wishart+Bartlett} reaches nominal — the residual ~10–15% Type I comes from the flexibility gap.

3. **Three-corrections picture** is documented in:
   - `claude/notes/2026-04-26-chisquare-recalibration.md` — earlier (variance-inflation) framing, partially superseded
   - `~/.claude/projects/.../memory/chisquare_three_corrections.md` — current understanding
   - `claude/derivations/mi_deviance_bias_derivation_v4.qmd` Sections 10 and 11 — formal statement

### Numerical target

At $N=250$, $q=9$, $p_{\text{analysis}}=32$, $p_{\text{sat}}=54$ ($\Delta p = 22$), mr=0.40, $M=50$:

$$\Delta(p_{\text{imp}}=54, p_{\text{analysis}}=32, \text{mr}=0.40, N=250) \;\approx\; +5 \text{ chi-square units}$$

This is what your derivation should reproduce in this regime. Cross-validate against:

- the saturated-imputation cells in `hpc/results-amelia-empri0.1/` (different N values: 100, 500, 1000)
- the PMM cells in `hpc/results-decomp/` (N: 100, 250, 1000) and `hpc/results-full-M100/` (full grid at $M=100$)

For each, the empirical $\Delta$ can be read off as $(\mathbb{E}[\chi^2_{\text{adhoc}}] - \mathbb{E}[\chi^2_{\text{com}}])$ minus the Bartlett residual.

---

## Mathematical framing

### Setup

- Observed-variable space: $\mathbb{R}^q$, $q$ continuous indicators
- Population: $Y \sim N(\mu, \Sigma)$ where $\Sigma = \Lambda \Lambda^T + \Psi$ has structured form (e.g., $M_1$ CFA)
- Missing pattern: MCAR with rate mr
- Imputation model: parametric family $\{p(Y \mid \phi) : \phi \in \mathbb{R}^{p_{\text{imp}}}\}$ with $\hat\phi$ FIML MLE on amputed data
- Analysis model: $\{p(Y \mid \theta) : \theta \in \mathbb{R}^{p_{\text{analysis}}}\}$ where $\Theta_{\text{analysis}} \subset \Phi_{\text{imp}}$ (analysis nests in imputation; this is what saturated imputation guarantees)

### Quantities

- $\chi^2_{\text{com}}(M)$: ML chi-square of $M$ vs $M_{\text{sat}}$ on complete data — well-studied
- $\chi^2_{\text{adhoc}}(M)$: average per-imputation chi-square, $\bar{D}(M) - \bar{D}(M_{\text{sat}})$ — what users compute by default
- $\chi^2_{\text{MI}}(M)$: bias-corrected via $\text{tr}(\text{RIV})$ difference — what the v4 framework prescribes

### What's known under congenial imputation ($\Theta_{\text{analysis}} = \Phi_{\text{imp}}$)

By v4's leading-order theory + Wishart correction:
$$\mathbb{E}[\chi^2_{\text{MI}, \text{Wishart}}] = \mathbb{E}[\chi^2_{\text{com}}] + O(N^{-2})$$

### What we need to derive: $\Delta(\cdot)$ under saturated imputation ($\Phi_{\text{imp}} \supsetneq \Theta_{\text{analysis}}$)

We want the closed form for:
$$\Delta = \mathbb{E}[\chi^2_{\text{adhoc, sat-impute}}] - \mathbb{E}[\chi^2_{\text{com}}] - \text{(Bartlett correction)} - \text{(Wishart contribution from tr(RIV) difference)}$$

### Conceptual approach (initial guess — verify or refute)

The intuition: when imputation has more parameters than analysis, each imputed dataset has variance in directions that the analysis model doesn't explain. Fitting $M_{\text{sat}}$ picks up some of this; fitting $M$ misses it; the difference inflates $\chi^2(M, M_{\text{sat}})$.

This sounds like a **non-centrality** problem: the asymptotic distribution of $\chi^2_{\text{adhoc}}$ under $H_0$ becomes non-central chi-square with non-centrality parameter $\lambda(\Phi_{\text{imp}} \setminus \Theta_{\text{analysis}}, \text{mr}, N)$ rather than central chi-square.

A natural starting point: characterize $\Delta$ as

$$\Delta(\cdot) \;\approx\; (p_{\text{imp}} - p_{\text{analysis}}) \cdot g(\text{mr}, N)$$

where $g$ is some function of the missing-data information loss. Empirically at our point: $\Delta = 5$, $p_{\text{imp}} - p_{\text{analysis}} = 22$, so $g \approx 0.23$. Whether $g$ is $\text{mr}/(2(1-\text{mr}))$ or something more complex is what the derivation needs to nail down.

### Plausible mathematical machinery

- **Non-central chi-square** asymptotic theory (Cudeck, Browne, Yuan & Bentler for SEM)
- **Misspecification literature** (Satorra-Saris, MacCallum 1990, Curran et al. 1996) — analyzes inflation when the analysis model is misspecified, but the *structure* of the misspecification under MI is different (the imputed data carries the imputation model's full $\phi$, not just $\theta$)
- **Trace identities** in multivariate Gaussian information matrices, $\mathcal{I}_{\text{imp}} = \mathcal{I}_{\text{analysis}} \oplus \mathcal{I}_{\Phi \setminus \Theta}$
- **Yuan, Marshall, Bentler (2003)** "On the asymptotic distribution of the test statistic when the data are missing at random" — closely related
- **Reilly (1995)** "Data analysis using hot deck multiple imputation" — discussion of chi-square inflation under non-congenial imputation

### Asymptotic regime to consider first

To make progress, start with the **clean asymptotic case**:
- $N \to \infty$, $M \to \infty$
- $q$ fixed, structured analysis model with $p_{\text{analysis}}$ params
- Saturated imputation: $\hat\phi$ is the unrestricted FIML estimator with $p_{\text{imp}} = q + q(q+1)/2$
- MCAR with fixed mr

The leading-order $\Delta$ in this regime is the cleanest target.

---

## Files to read (in priority order)

1. **`claude/derivations/mi_deviance_bias_derivation_v4.qmd`** — current state of theory, especially Sections 10 and 11. Read for what assumptions v4 makes and where they break.
2. **`claude/notes/2026-04-26-chisquare-recalibration.md`** — earlier write-up of the empirical findings; some framings have been refined (see chisquare_three_corrections.md memory).
3. **`hpc/check_wishart_chi2.R`** + **`hpc/check_bartlett_wishart.R`** — exactly how Wishart and Bartlett corrections were applied empirically. Useful for matching empirical $\Delta$ values.
4. **`hpc/check_mvn_M1_chi2.R`** — true-model imputation comparison, which is the "no flexibility gap" baseline.
5. **`packages/miicsem/R/run_replication.R`** — how `chi2_MI` is actually computed in the simulation pipeline.

## Empirical artifacts to validate against

Use these for cross-cell verification of your derived $\Delta(\cdot)$:

- `hpc/results-decomp/results_combined.rds` — PMM, M=50, N ∈ {100, 250, 1000}, mr=0.40
- `hpc/results-amelia-empri0.1/results_combined.rds` — saturated MVN imputation (proper, mild empri), M=100, N ∈ {100, 500, 1000}, mr ∈ {0.20, 0.40}
- `hpc/results-amelia-empri100/results_combined.rds` — same but with empri=100*N (severely uncongenial)
- `hpc/results-mvn_M1-N250/` and `hpc/results-mvn_M1-N1000/` — true-model imputation, congenial baseline
- `hpc/results-full-M100/results_combined.rds` — full PMM grid, $M=100$, the original Study 2 baseline (no decomp columns but has $\chi^2$ values)

For each cell, the relevant empirical observable is $\mathbb{E}[\chi^2_{\text{adhoc}}] - \text{df}$ minus a Bartlett residual.

## Deliverables

1. **`claude/derivations/flexibility_gap_derivation_v0.qmd`** — first-pass analytic derivation of $\Delta(p_{\text{imp}}, p_{\text{analysis}}, \text{mr}, N)$
2. Cross-validation table: predicted $\Delta$ from your formula vs empirical $\Delta$ from each available cell
3. If the closed form is clean, propose appending it as **Section 11.5: Flexibility-Gap Correction** to v4.qmd
4. Update `chisquare_three_corrections.md` memory and add a manuscript-ready statement of the corrected $\chi^2$ formula incorporating all three terms

## Acceptance criteria

- Derived $\Delta$ matches empirical $\Delta$ within 1 chi-square unit (Monte Carlo SE) at the **mvn_M1** cells (where $\Delta = 0$ is expected) AND at least one saturated-imputation cell.
- Functional form has a sensible limit: $\Delta \to 0$ as $p_{\text{imp}} \to p_{\text{analysis}}$, and $\Delta$ scales monotonically with mr.
- Derivation is rigorous to leading order in $N$ (matches Bartlett's level of rigor), with $O(N^{-1})$ residual stated.

## Useful identities and notational conventions

Use the v4 conventions throughout:

- $\theta$ for analysis-model parameters, $\phi$ (or $\tilde\phi$ for imputation draws) for imputation-model parameters
- $\bar Q_{MI}$ for Q-function evaluated at $\bar\theta_{MI}$ (Rubin's rule pooled mean)
- $\text{RIV}_M = (1+1/M) W_M^{-1} B_M$ for the relative-increase-in-variance matrix
- $\hat V_{\text{obs},M}$ and $\hat V_{\text{com},M}$ for FIML and complete-data inverse-Hessian variance estimates
- Wishart factor: $w(n, p) = (n - p - 1)/n$
- Bartlett factor: $b(N, p) \approx 1/(1 + c(p)/N)$

## Important caveat

The two earlier derivation revisions (v3, v4) had subtle errors that took multiple iterations to surface. Be skeptical of any "clean closed form" until verified against the empirical cells. Plan to commit a v0 first pass and let it be wrong before committing a v1.

## Coordinate with the codebase

- The `miicsem` R package (version 0.5.2 at handoff) currently exposes raw `tr_RIV`. If your derivation needs new sample-statistic outputs, propose them as an extension; do not modify the existing simulation results.
- All sims used master_seed=32897891 with stage-split SeedMaker per rep. Reproduction is one-command via `miicsem::run_simulation()`.
