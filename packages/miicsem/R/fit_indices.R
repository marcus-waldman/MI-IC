#' Compute MI-Corrected SEM Fit Indices (CFI, TLI, RMSEA)
#'
#' For each candidate analysis model (rows of \code{chi2_df} other than
#' "Mnull" and "Msat"), computes three variants of CFI, TLI, and RMSEA:
#' the oracle complete-data variant (\code{_com}), the standard ad-hoc
#' MI-pooled variant (\code{_adhoc}), and the v4.3 MI-corrected variant
#' (\code{_MI}). The MI variant uses the \code{tr(RIV)}-corrected
#' \code{chi2_MI} statistic, which is mean-calibrated to \code{chi2_com}
#' under saturated proper imputation
#' (see \code{claude/derivations/mi_deviance_bias_derivation_v4.qmd}
#' Sections 9 and 10).
#'
#' Index formulas (Bentler 1990; Tucker & Lewis 1973; Steiger & Lind 1980):
#'
#' \itemize{
#'   \item \strong{CFI:}
#'     \eqn{1 - \max(\chi^2(M) - df_M, 0) /
#'          \max(\chi^2(M) - df_M,\, \chi^2(M_0) - df_0,\, 0)}
#'   \item \strong{TLI} (a.k.a.\ NNFI):
#'     \eqn{[\chi^2(M_0)/df_0 - \chi^2(M)/df_M] / [\chi^2(M_0)/df_0 - 1]}
#'   \item \strong{RMSEA:}
#'     \eqn{\sqrt{\max(\chi^2(M) - df_M, 0) / [df_M (N - 1)]}}
#' }
#'
#' Each is computed in three variants by substituting \code{chi2_com},
#' \code{chi2_adhoc}, or \code{chi2_MI} for \eqn{\chi^2(\cdot)} (the same
#' substitution applies to both candidate and null in CFI/TLI).
#'
#' @param chi2_df data.frame from \code{\link{compute_chi_squares}}; rows
#'   are model names with columns \code{chi2_com}, \code{chi2_adhoc},
#'   \code{chi2_MI}, \code{df}. Must include an "Mnull" row when CFI/TLI
#'   are requested; otherwise CFI and TLI columns will be \code{NA}.
#' @param n Sample size (numeric). Used by RMSEA.
#' @return data.frame with one row per candidate model (rows of
#'   \code{chi2_df} other than "Mnull" and "Msat"; "Msat" is not in
#'   \code{chi2_df} by construction). Columns:
#'   \code{cfi_com}, \code{cfi_adhoc}, \code{cfi_MI},
#'   \code{tli_com}, \code{tli_adhoc}, \code{tli_MI},
#'   \code{rmsea_com}, \code{rmsea_adhoc}, \code{rmsea_MI}.
#' @export
compute_fit_indices <- function(chi2_df, n) {

  null_name <- "Mnull"
  has_null <- null_name %in% rownames(chi2_df)

  if (has_null) {
    chi_null    <- chi2_df[null_name, ]
    chi0_com    <- chi_null$chi2_com
    chi0_adhoc  <- chi_null$chi2_adhoc
    chi0_MI     <- chi_null$chi2_MI
    df0         <- chi_null$df
  } else {
    chi0_com <- chi0_adhoc <- chi0_MI <- df0 <- NA_real_
  }

  candidates <- setdiff(rownames(chi2_df), null_name)

  cfi_one <- function(chi_M, df_M, chi_0, df_0) {
    if (is.na(chi_M) || is.na(df_M) || is.na(chi_0) || is.na(df_0))
      return(NA_real_)
    num <- max(chi_M - df_M, 0)
    den <- max(chi_M - df_M, chi_0 - df_0, 0)
    if (den <= 0) return(1.0)
    1 - num / den
  }

  tli_one <- function(chi_M, df_M, chi_0, df_0) {
    if (is.na(chi_M) || is.na(df_M) || is.na(chi_0) || is.na(df_0) ||
        df_M <= 0 || df_0 <= 0)
      return(NA_real_)
    num <- chi_0 / df_0 - chi_M / df_M
    den <- chi_0 / df_0 - 1
    if (abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
  }

  rmsea_one <- function(chi_M, df_M, n_val) {
    if (is.na(chi_M) || is.na(df_M) || df_M <= 0 || n_val <= 1)
      return(NA_real_)
    num <- max(chi_M - df_M, 0)
    sqrt(num / (df_M * (n_val - 1)))
  }

  rows <- lapply(candidates, function(mname) {
    r       <- chi2_df[mname, ]
    chi_com <- r$chi2_com;   chi_adhoc <- r$chi2_adhoc;   chi_MI <- r$chi2_MI
    df_M    <- r$df

    data.frame(
      cfi_com    = cfi_one(chi_com,   df_M, chi0_com,   df0),
      cfi_adhoc  = cfi_one(chi_adhoc, df_M, chi0_adhoc, df0),
      cfi_MI     = cfi_one(chi_MI,    df_M, chi0_MI,    df0),
      tli_com    = tli_one(chi_com,   df_M, chi0_com,   df0),
      tli_adhoc  = tli_one(chi_adhoc, df_M, chi0_adhoc, df0),
      tli_MI     = tli_one(chi_MI,    df_M, chi0_MI,    df0),
      rmsea_com  = rmsea_one(chi_com,   df_M, n),
      rmsea_adhoc = rmsea_one(chi_adhoc, df_M, n),
      rmsea_MI   = rmsea_one(chi_MI,    df_M, n),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- candidates
  out
}
