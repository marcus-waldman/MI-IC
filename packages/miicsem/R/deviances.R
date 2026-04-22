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
#' @param decomp_logliks Optional data.frame from
#'   \code{\link{compute_decomp_logliks}} with columns \code{model},
#'   \code{loglik_com_at_obs}, \code{loglik_com_at_pooled}. If supplied,
#'   these columns plus \code{loglik_com_at_com} and
#'   \code{mean_loglik_at_pooled} are appended to the output so the
#'   three-term bias decomposition can be computed offline.
#' @return data.frame with one row per model and columns
#'   \code{DEV_com}, \code{DEV_adhoc}, \code{MI_DEVIANCE},
#'   \code{MR_DEVIANCE} (NA here — populated later),
#'   \code{npar}, \code{tr_RIV}, \code{tr_RIV_fiml}, plus (optionally)
#'   \code{loglik_com_at_com}, \code{loglik_com_at_obs},
#'   \code{loglik_com_at_pooled}, \code{mean_loglik_at_pooled}.
compute_deviances <- function(complete_fits, mi_fits,
                              observed_fits = NULL, n,
                              decomp_logliks = NULL) {
  model_names <- names(complete_fits)

  rows <- lapply(model_names, function(mname) {
    cf <- complete_fits[[mname]]
    mf <- mi_fits[[mname]]
    of <- if (!is.null(observed_fits)) observed_fits[[mname]] else NULL

    DEV_com <- if (!is.na(cf$loglik) && isTRUE(cf$converged)) {
      -2 * cf$loglik
    } else NA_real_

    if (isTRUE(mf$success)) {
      DEV_adhoc    <- -2 * mf$mean_loglik
      logliks_pool <- mf$logliks_at_pooled
      logliks_pool <- logliks_pool[!is.na(logliks_pool)]
      mean_loglik_at_pooled <- if (length(logliks_pool) > 0) {
        mean(logliks_pool)
      } else NA_real_
      MI_DEVIANCE <- if (!is.na(mean_loglik_at_pooled)) {
        -2 * mean_loglik_at_pooled + mf$tr_RIV
      } else NA_real_
      npar   <- mf$npar
      tr_RIV <- mf$tr_RIV
    } else {
      DEV_adhoc             <- NA_real_
      MI_DEVIANCE           <- NA_real_
      mean_loglik_at_pooled <- NA_real_
      npar                  <- if (!is.na(cf$npar)) cf$npar else NA_integer_
      tr_RIV                <- NA_real_
    }

    tr_RIV_fiml <- NA_real_
    if (!is.null(of) && isTRUE(of$converged) &&
        !is.null(of$vcov) && !is.null(cf$vcov)) {
      V_obs <- of$vcov
      V_com <- cf$vcov
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

    # Decomposition log-likelihoods (oracle ell_com(theta_com) comes
    # straight from the complete-data fit; the other two come from
    # decomp_logliks if supplied).
    loglik_com_at_com    <- if (!is.na(cf$loglik) && isTRUE(cf$converged)) {
      as.numeric(cf$loglik)
    } else NA_real_
    loglik_com_at_obs    <- NA_real_
    loglik_com_at_pooled <- NA_real_
    if (!is.null(decomp_logliks)) {
      r <- decomp_logliks[decomp_logliks$model == mname, ]
      if (nrow(r) == 1) {
        loglik_com_at_obs    <- r$loglik_com_at_obs
        loglik_com_at_pooled <- r$loglik_com_at_pooled
      }
    }

    data.frame(
      DEV_com               = DEV_com,
      DEV_adhoc             = DEV_adhoc,
      MI_DEVIANCE           = MI_DEVIANCE,
      MR_DEVIANCE           = NA_real_,
      npar                  = npar,
      tr_RIV                = tr_RIV,
      tr_RIV_fiml           = tr_RIV_fiml,
      loglik_com_at_com     = loglik_com_at_com,
      loglik_com_at_obs     = loglik_com_at_obs,
      loglik_com_at_pooled  = loglik_com_at_pooled,
      mean_loglik_at_pooled = mean_loglik_at_pooled,
      stringsAsFactors      = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- model_names
  out
}
