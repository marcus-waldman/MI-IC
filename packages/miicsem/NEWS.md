# miicsem 0.5.3 (2026-04-28)

## New features

- **Null/independence model**: new `get_null_model(var_names)` returns
  the lavaan syntax for the independence baseline (free variances,
  fixed-zero covariances, free intercepts via `meanstructure = TRUE`).
  For the Bollen 9-indicator design, `p_M0 = 18` and `df_0 = 36` against
  the saturated reference.
- **MI-corrected SEM fit indices**: new `compute_fit_indices(chi2_df, n)`
  returns CFI, TLI, and RMSEA in three variants (oracle `_com`, ad-hoc
  `_adhoc`, and v4.3 MI-corrected `_MI`). The MI variant is mean-calibrated
  to the oracle under saturated proper imputation; see
  `claude/derivations/mi_deviance_bias_derivation_v4.qmd` Sections 9-10.
- `run_one_rep()` now fits the null model alongside the saturated reference
  on every imputation and returns a new `fit_indices_df` element in its
  output list.

## Behavior changes

- `compute_chi_squares()` now adds a `chi2_adhoc` column to the returned
  `chi2_df` (previously only `chi2_com`, `chi2_MI`, `chi2_D3` were stored).
  This is needed by `compute_fit_indices()` and is generally useful for
  downstream analysis of the bias.
- The `chi2_df` output of `run_one_rep()` now includes a row for
  `Mnull` (in addition to the candidate models). The IC pipeline
  (`compute_all_models_ic`) is unaffected — it operates only on the
  candidate model list, not on the chi-square baselines.

## Compatibility

- Existing simulation result `.rds` files from miicsem ≤ 0.5.2 do not
  contain `fit_indices_df` and will work unchanged with downstream
  analysis scripts that don't reference it.
- Runs initiated with miicsem 0.5.3 will include `fit_indices_df` in
  every per-rep output and an extra `Mnull` row in `chi2_df` and `dev_df`.

# miicsem 0.5.2 (2026-04-26)

- Added `mvn_M1` imputation method (truly congenial: fits M1 by FIML,
  draws Y_mis from M1's implied conditional). Used to verify the v4
  Term 1 result under the strictest possible congenial regime.

# miicsem 0.5.1 (2026-04-12)

- Moved Amelia from Suggests to Imports.

# miicsem 0.5.0 (2026-04-08)

- Amelia × empri sweep added: `mice_method = "amelia"` with
  `amelia_empri_frac` parameter for the congeniality dose-response
  simulation.
