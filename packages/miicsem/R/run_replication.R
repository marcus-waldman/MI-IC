#' Run One Replication of the SEM Simulation
#'
#' Generates complete data, fits all 12 oracle models, amputes under MCAR,
#' runs MI via \code{mice}, fits and pools all 12 models, computes the
#' seven IC values, and returns selections and diagnostics.
#'
#' @param rep_id Integer replication ID (for logging).
#' @param n Sample size.
#' @param miss_rate Proportion of incomplete cases.
#' @param config Configuration list from \code{\link{get_config}}.
#' @param seed_data Seed for data generation.
#' @param seed_ampute Seed for amputation.
#' @param seed_impute Seed for mice imputation.
#' @param pop_starts Optional named list of warm-start coefficient
#'   vectors covering all entries in \code{models_with_sat} (see
#'   \code{\link{compute_pop_starts}}).
#' @return List with \code{ic_df}, \code{selections},
#'   \code{convergence}, \code{tr_RIVs}; or NULL on imputation failure.
#' @export
run_one_rep <- function(rep_id, n, miss_rate, config,
                        seed_data, seed_ampute, seed_impute,
                        pop_starts = NULL) {

  models <- get_sim1_models()
  models_with_sat <- c(models, list(Msat = get_saturated_model(config$var_names)))

  data_complete <- generate_complete_data(n, config, seed_data)

  complete_fits <- fit_complete(data_complete, models_with_sat, pop_starts = pop_starts)

  data_miss <- ampute_data(data_complete, miss_rate, seed_ampute)

  # FIML fits on amputed data → V_obs for tr(RIV_fiml) reference
  observed_fits <- fit_observed(data_miss, models_with_sat,
                                pop_starts = pop_starts)

  # Dispatch on imputation method.  "pmm" and "norm" go through mice;
  # "mvn_msat" uses the congenial saturated-MVN imputer that matches the
  # analysis model's multivariate normal family.
  mice_method <- config$mice_method
  if (is.null(mice_method)) mice_method <- "pmm"

  if (mice_method %in% c("pmm", "norm", "norm.predict", "norm.nob",
                         "norm.boot")) {
    imp <- tryCatch({
      mice::mice(
        data_miss,
        m         = config$M,
        method    = mice_method,
        maxit     = config$mice_maxit,
        seed      = seed_impute,
        printFlag = FALSE
      )
    }, error = function(e) NULL)

    if (is.null(imp)) return(NULL)

    imputed_list <- lapply(seq_len(config$M), function(m) {
      mice::complete(imp, action = m)
    })
  } else if (mice_method == "mvn_msat") {
    imputed_list <- tryCatch({
      impute_from_saturated_mvn(
        data_miss = data_miss,
        M         = config$M,
        seed      = seed_impute,
        var_names = config$var_names
      )
    }, error = function(e) NULL)
    if (is.null(imputed_list)) return(NULL)
  } else if (mice_method == "amelia") {
    # Amelia's EMB: bootstrap the data, EM to get (mu*, Sigma*),
    # draw Y_mis | Y_obs from the resulting joint MVN.  Proper MI with
    # joint-coherent imputation model.
    set.seed(seed_impute)
    imp_obj <- tryCatch({
      utils::capture.output({
        a <- Amelia::amelia(
          x     = data_miss,
          m     = config$M,
          p2s   = 0,
          boot.type = "ordinary"
        )
      }, type = "message")
      a
    }, error = function(e) NULL)
    if (is.null(imp_obj) || inherits(imp_obj, "try-error")) return(NULL)
    # amelia returns NULLs in $imputations for any failed imputation
    imputed_list <- imp_obj$imputations
    imputed_list <- imputed_list[!vapply(imputed_list, is.null, logical(1))]
    if (length(imputed_list) < 2) return(NULL)
  } else {
    stop(sprintf("Unknown mice_method: '%s'", mice_method))
  }

  mi_fits <- fit_mi_models(imputed_list, models_with_sat, config, pop_starts = pop_starts)

  # Three-term bias decomposition: evaluate complete-data log-likelihood
  # at theta_obs_fiml and at theta_pooled.  Together with the oracle
  # ell_com(theta_com) (already in complete_fits) and Q_bar(theta_pooled)
  # (in mi_fits$logliks_at_pooled), this lets us decompose the pooled-MI
  # vs oracle bias into three per-replicate pieces:
  #   Term 1 (Imputation) = Q_bar(theta_pooled) - ell_com(theta_pooled)
  #   Term 2 (Pooling)    = ell_com(theta_pooled) - ell_com(theta_obs)
  #   Term 3 (Estimation) = ell_com(theta_obs)    - ell_com(theta_com)
  decomp_logliks <- compute_decomp_logliks(
    complete_fits  = complete_fits,
    observed_fits  = observed_fits,
    mi_fits        = mi_fits,
    data_complete  = data_complete
  )

  # Deviances (all on -2 log-likelihood scale) for every model incl. Msat
  dev_df <- compute_deviances(complete_fits, mi_fits, observed_fits, n,
                              decomp_logliks = decomp_logliks)

  # Chi-squares (candidate models only) vs saturated reference; fills in
  # MR_DEVIANCE column of dev_df.
  chi2_res <- compute_chi_squares(complete_fits, mi_fits, dev_df)
  chi2_df  <- chi2_res$chi2_df
  dev_df   <- chi2_res$dev_df

  # Seven IC methods for candidate models only
  ic_df <- compute_all_models_ic(complete_fits[names(models)],
                                 mi_fits[names(models)], n)

  selections <- select_models(ic_df)

  convergence <- vapply(mi_fits, function(f) {
    if (!f$success) return(0L)
    as.integer(f$converged_count)
  }, integer(1))

  tr_RIVs <- vapply(mi_fits, function(f) {
    if (!f$success) return(NA_real_)
    f$tr_RIV
  }, numeric(1))

  list(
    rep_id      = rep_id,
    n           = n,
    miss_rate   = miss_rate,
    ic_df       = ic_df,
    dev_df      = dev_df,
    chi2_df     = chi2_df,
    selections  = selections,
    convergence = convergence,
    tr_RIVs     = tr_RIVs
  )
}
