# ============================================================================
# Single Replication Wrapper for Study 2 (SEM)
# ============================================================================

#' Run One Replication of the SEM Simulation
#'
#' @param rep_id Integer replication ID (for logging)
#' @param n Sample size
#' @param miss_rate Proportion of incomplete cases
#' @param config Configuration list from get_config()
#' @param seed_data Seed for data generation
#' @param seed_ampute Seed for amputation
#' @param seed_impute Seed for mice imputation
#' @return List with: ic_df, selections, diagnostics; or NULL on failure
run_one_rep <- function(rep_id, n, miss_rate, config,
                        seed_data, seed_ampute, seed_impute) {

  models <- get_sim1_models()

  # 1. Generate complete data
  data_complete <- generate_complete_data(n, config, seed_data)

  # 2. Fit complete-data models (oracle)
  complete_fits <- fit_complete(data_complete, models)

  # 3. Ampute
  data_miss <- ampute_data(data_complete, miss_rate, seed_ampute)

  # 4. Multiple imputation with mice
  imp <- tryCatch({
    mice::mice(
      data_miss,
      m       = config$M,
      method  = config$mice_method,
      maxit   = config$mice_maxit,
      seed    = seed_impute,
      printFlag = FALSE
    )
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(imp)) return(NULL)

  # Extract imputed datasets as a list
  imputed_list <- lapply(seq_len(config$M), function(m) {
    mice::complete(imp, action = m)
  })

  # 5. Fit all models to MI data and pool
  mi_fits <- fit_mi_models(imputed_list, models, config)

  # 6. Compute IC for all models
  ic_df <- compute_all_models_ic(complete_fits, mi_fits, n)

  # 7. Select best model per IC method
  selections <- select_models(ic_df)

  # 8. Diagnostics
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
