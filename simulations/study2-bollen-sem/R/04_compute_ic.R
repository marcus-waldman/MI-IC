# ============================================================================
# Information Criteria Computation for Study 2 (SEM)
# ============================================================================
# Five IC methods (AIC + BIC variants), D_LR-based AIC deferred.
# ============================================================================

#' Compute All IC Methods for a Single Model
#'
#' @param complete_fit List with loglik, npar from fit_complete()
#' @param mi_fit List with mean_loglik, npar, tr_RIV, logliks_at_pooled from fit_mi_models()
#' @param n Sample size
#' @return Named numeric vector of IC values (NA if inputs are invalid)
compute_all_ic <- function(complete_fit, mi_fit, n) {

  out <- c(
    AIC_com   = NA_real_, BIC_com   = NA_real_,
    AIC_adhoc = NA_real_, BIC_adhoc = NA_real_,
    AICcd     = NA_real_,
    MI_AIC    = NA_real_, MI_BIC    = NA_real_
  )

  # --- 1. Complete-data IC (oracle) ---
  if (!is.na(complete_fit$loglik) && complete_fit$converged) {
    p <- complete_fit$npar
    out["AIC_com"] <- -2 * complete_fit$loglik + 2 * p
    out["BIC_com"] <- -2 * complete_fit$loglik + p * log(n)
  }

  # --- MI-based ICs require successful MI fit ---
  if (!mi_fit$success) return(out)

  p      <- mi_fit$npar
  tr_RIV <- mi_fit$tr_RIV

  # --- 2. Ad hoc IC: mean of per-imputation IC ---
  # mean_loglik is mean of logliks at each imputation's own MLE
  out["AIC_adhoc"] <- -2 * mi_fit$mean_loglik + 2 * p
  out["BIC_adhoc"] <- -2 * mi_fit$mean_loglik + p * log(n)

  # --- 3. AICcd (Cavanaugh & Shumway): evaluated at own MLEs, penalty 2p + 2*tr(RIV) ---
  out["AICcd"] <- -2 * mi_fit$mean_loglik + 2 * p + 2 * tr_RIV

  # --- 4. MI-AIC and MI-BIC (our proposal): evaluated at pooled theta ---
  logliks_pooled <- mi_fit$logliks_at_pooled
  logliks_pooled <- logliks_pooled[!is.na(logliks_pooled)]
  if (length(logliks_pooled) > 0) {
    mean_loglik_pooled <- mean(logliks_pooled)
    out["MI_AIC"] <- -2 * mean_loglik_pooled + 2 * p + tr_RIV
    out["MI_BIC"] <- -2 * mean_loglik_pooled + p * log(n) + tr_RIV
  }

  return(out)
}


#' Compute IC for All 12 Models
#'
#' @param complete_fits Named list from fit_complete()
#' @param mi_fits Named list from fit_mi_models()
#' @param n Sample size
#' @return data.frame with model names as rows and IC values as columns
compute_all_models_ic <- function(complete_fits, mi_fits, n) {
  model_names <- names(complete_fits)

  ic_mat <- t(sapply(model_names, function(mname) {
    compute_all_ic(complete_fits[[mname]], mi_fits[[mname]], n)
  }))

  as.data.frame(ic_mat)
}


#' Select Best Model per IC Method
#'
#' @param ic_df data.frame from compute_all_models_ic()
#' @return Named character vector: method -> selected model name
select_models <- function(ic_df) {
  vapply(colnames(ic_df), function(method) {
    vals <- ic_df[[method]]
    if (all(is.na(vals))) return(NA_character_)
    rownames(ic_df)[which.min(vals)]
  }, character(1))
}
