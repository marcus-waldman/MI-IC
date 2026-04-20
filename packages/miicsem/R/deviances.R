#' Compute Per-Model Log-Likelihood / Deviance Quantities
#'
#' All four quantities are on the deviance (-2 log-likelihood) scale.
#'
#' \itemize{
#'   \item \code{DEV_com}: oracle; \eqn{-2\ell_{\text{com}}(\hat\theta_{\text{com}})}
#'         from the complete-data fit.
#'   \item \code{DEV_adhoc}: ad hoc; mean over imputations of
#'         \eqn{-2\ell_m(\hat\theta_m)} (each imputation at its own MLE).
#'   \item \code{MI_DEVIANCE}: our debiased pooled deviance,
#'         \eqn{-2\bar Q_{M_j}(\bar\theta_{M_j}) + \text{tr}(RIV_{M_j})}.
#'   \item \code{MR_DEVIANCE}: Meng-Rubin (1992) anchored at saturated;
#'         computed separately by \code{\link{compute_chi_squares}} since
#'         it depends on both the candidate and saturated fits.
#' }
#'
#' This function populates \code{DEV_com}, \code{DEV_adhoc}, and
#' \code{MI_DEVIANCE} for every model (candidates + saturated).
#' \code{MR_DEVIANCE} is filled in later by
#' \code{\link{compute_chi_squares}}.
#'
#' @param complete_fits Named list from \code{\link{fit_complete}}.
#' @param mi_fits Named list from \code{\link{fit_mi_models}}.
#' @param n Sample size (unused here but kept for signature consistency).
#' @return data.frame with one row per model and columns
#'   \code{DEV_com}, \code{DEV_adhoc}, \code{MI_DEVIANCE},
#'   \code{MR_DEVIANCE} (initialised to NA — populated later),
#'   \code{npar}, \code{tr_RIV}.
compute_deviances <- function(complete_fits, mi_fits, n) {
  model_names <- names(complete_fits)

  rows <- lapply(model_names, function(mname) {
    cf <- complete_fits[[mname]]
    mf <- mi_fits[[mname]]

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
      DEV_adhoc <- NA_real_
      MI_DEVIANCE <- NA_real_
      npar <- if (!is.na(cf$npar)) cf$npar else NA_integer_
      tr_RIV <- NA_real_
    }

    data.frame(
      DEV_com     = DEV_com,
      DEV_adhoc   = DEV_adhoc,
      MI_DEVIANCE = MI_DEVIANCE,
      MR_DEVIANCE = NA_real_,
      npar        = npar,
      tr_RIV      = tr_RIV,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- model_names
  out
}
