#' Compute Per-Model Log-Likelihood / Deviance / RIV Quantities
#'
#' All deviances are on the -2 log-likelihood scale.
#'
#' \itemize{
#'   \item \code{DEV_com}: oracle; \eqn{-2\ell_{\text{com}}(\hat\theta_{\text{com}})}
#'         from the complete-data fit.
#'   \item \code{DEV_adhoc}: ad hoc; mean over imputations of
#'         \eqn{-2\ell_m(\hat\theta_m)} (each imputation at its own MLE).
#'   \item \code{MI_DEVIANCE}: our debiased pooled deviance,
#'         \eqn{-2\bar Q_{M_j}(\bar\theta_{M_j}) + \text{tr}(RIV_{M_j})}.
#'   \item \code{MR_DEVIANCE}: Meng-Rubin (1992) anchored at saturated;
#'         filled in later by \code{\link{compute_chi_squares}}.
#'   \item \code{tr_RIV}: MI-based (1 + 1/M) tr(W^{-1} B).
#'   \item \code{tr_RIV_fiml}: theoretical asymptotic target,
#'         \eqn{\operatorname{tr}(\hat V_{\text{obs}} \hat V_{\text{com}}^{-1}) - p},
#'         evaluated at each fit's own MLE (V_com at theta_hat_com,
#'         V_obs at theta_hat_fiml). Consistent for the asymptotic RIV
#'         under MAR.
#' }
#'
#' @param complete_fits Named list from \code{\link{fit_complete}}; each
#'   element must contain \code{vcov} (the complete-data vcov).
#' @param mi_fits Named list from \code{\link{fit_mi_models}}.
#' @param observed_fits Named list from \code{\link{fit_observed}}; each
#'   element must contain \code{vcov} (FIML vcov). Optional — if NULL,
#'   the \code{tr_RIV_fiml} column is filled with NA.
#' @param n Sample size (unused here; kept for signature symmetry).
#' @return data.frame with one row per model and columns
#'   \code{DEV_com}, \code{DEV_adhoc}, \code{MI_DEVIANCE},
#'   \code{MR_DEVIANCE} (NA here — populated later),
#'   \code{npar}, \code{tr_RIV}, \code{tr_RIV_fiml}.
compute_deviances <- function(complete_fits, mi_fits,
                              observed_fits = NULL, n) {
  model_names <- names(complete_fits)

  rows <- lapply(model_names, function(mname) {
    cf <- complete_fits[[mname]]
    mf <- mi_fits[[mname]]
    of <- if (!is.null(observed_fits)) observed_fits[[mname]] else NULL

    DEV_com <- if (!is.na(cf$loglik) && isTRUE(cf$converged)) {
      -2 * cf$loglik
    } else NA_real_

    if (isTRUE(mf$success)) {
      DEV_adhoc   <- -2 * mf$mean_loglik
      logliks_pool <- mf$logliks_at_pooled
      logliks_pool <- logliks_pool[!is.na(logliks_pool)]
      MI_DEVIANCE <- if (length(logliks_pool) > 0) {
        -2 * mean(logliks_pool) + mf$tr_RIV
      } else NA_real_
      npar   <- mf$npar
      tr_RIV <- mf$tr_RIV
    } else {
      DEV_adhoc   <- NA_real_
      MI_DEVIANCE <- NA_real_
      npar   <- if (!is.na(cf$npar)) cf$npar else NA_integer_
      tr_RIV <- NA_real_
    }

    tr_RIV_fiml <- NA_real_
    if (!is.null(of) && isTRUE(of$converged) &&
        !is.null(of$vcov) && !is.null(cf$vcov)) {
      V_obs <- of$vcov
      V_com <- cf$vcov
      # Match parameter order by name (FIML fit may have different order)
      common_names <- intersect(rownames(V_com), rownames(V_obs))
      if (length(common_names) >= 2 &&
          length(common_names) == nrow(V_com)) {
        V_com_sub <- V_com[common_names, common_names, drop = FALSE]
        V_obs_sub <- V_obs[common_names, common_names, drop = FALSE]
        p_sub <- length(common_names)
        tr_RIV_fiml <- tryCatch(
          sum(diag(V_obs_sub %*% solve(V_com_sub))) - p_sub,
          error = function(e) NA_real_
        )
      }
    }

    data.frame(
      DEV_com     = DEV_com,
      DEV_adhoc   = DEV_adhoc,
      MI_DEVIANCE = MI_DEVIANCE,
      MR_DEVIANCE = NA_real_,
      npar        = npar,
      tr_RIV      = tr_RIV,
      tr_RIV_fiml = tr_RIV_fiml,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- model_names
  out
}
