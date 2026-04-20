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
#' @return List with \code{ic_df}, \code{selections},
#'   \code{convergence}, \code{tr_RIVs}; or NULL on imputation failure.
#' @export
run_one_rep <- function(rep_id, n, miss_rate, config,
                        seed_data, seed_ampute, seed_impute) {

  models <- get_sim1_models()

  data_complete <- generate_complete_data(n, config, seed_data)

  complete_fits <- fit_complete(data_complete, models)

  data_miss <- ampute_data(data_complete, miss_rate, seed_ampute)

  imp <- tryCatch({
    mice::mice(
      data_miss,
      m         = config$M,
      method    = config$mice_method,
      maxit     = config$mice_maxit,
      seed      = seed_impute,
      printFlag = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(imp)) return(NULL)

  imputed_list <- lapply(seq_len(config$M), function(m) {
    mice::complete(imp, action = m)
  })

  mi_fits <- fit_mi_models(imputed_list, models, config)

  ic_df <- compute_all_models_ic(complete_fits, mi_fits, n)

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
    selections  = selections,
    convergence = convergence,
    tr_RIVs     = tr_RIVs
  )
}
