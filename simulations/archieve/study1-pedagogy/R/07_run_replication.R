# ============================================================================
# Run Single Replication for Study 1
# ============================================================================
# Wrapper function that runs the complete pipeline for one replication
# ============================================================================

#' Run Single Replication
#'
#' Executes the full simulation pipeline for one replication of one condition.
#' Automatically chooses between analytical (M = Inf) and Monte Carlo (M = finite)
#' approaches based on config$M.
#'
#' @param rep_id Integer. Replication number (1 to n_reps)
#' @param mechanism Character. "MCAR", "MAR", or "MNAR"
#' @param imputation_model Character. "matched", "saturated", or "uncongenial"
#' @param config List. Configuration from config.R. If config$M = Inf, uses
#'   analytical approach; otherwise uses brms Monte Carlo.
#' @param seeds List. Pre-generated seeds from setup_seeds()
#' @param verbose Logical. Print progress messages?
#' @return List with metrics from compute_study1_metrics(), or NULL if failed
run_single_replication <- function(rep_id,
                                    mechanism,
                                    imputation_model,
                                    config,
                                    seeds,
                                    verbose = FALSE) {

  # Get seeds for this replication and condition
  rep_seeds <- get_replication_seeds(seeds, rep_id, mechanism, imputation_model, config)

  # Get population parameters for comparison
  pop_params <- get_population_params(config)

  # Check if using analytical approach (M = Inf)
  use_analytical <- is.infinite(config$M)

  if (verbose) {
    method_str <- if (use_analytical) "ANALYTICAL" else sprintf("brms (M=%d)", config$M)
    cat(sprintf("Rep %d | %s | %s | %s | seed_data=%d\n",
                rep_id, mechanism, imputation_model, method_str, rep_seeds$seed_data))
  }

  result <- tryCatch({

    # 1. Generate complete data
    data_complete <- generate_mediation_data_from_config(config, seed = rep_seeds$seed_data)

    # 2. Ampute data
    ampute_result <- ampute_mediation(
      data = data_complete,
      prop_missing = config$prop_missing,
      mechanism = mechanism,
      seed = rep_seeds$seed_ampute
    )
    data_miss <- ampute_result$data_miss

    # Branch based on M = Inf vs finite
    if (use_analytical) {
      # ========== ANALYTICAL APPROACH (M = Inf) ==========
      # Note: For analytical approach, imputation_model is currently ignored

      # (assumes matched/congenial model). Future: implement uncongenial analytical.

      metrics <- compute_study1_metrics_analytical(
        data_complete = data_complete,
        data_miss = data_miss,
        pop_params = pop_params
      )

      # Add condition identifiers
      metrics$rep_id <- rep_id
      metrics$mechanism <- mechanism
      metrics$imputation_model <- imputation_model

      # No brms convergence info for analytical
      metrics$brms_converged <- NA
      metrics$brms_max_rhat <- NA

    } else {
      # ========== MONTE CARLO APPROACH (M = finite) ==========

      # 3. Fit brms imputation model
      brms_fit <- fit_brms_imputation_from_config(
        data_miss = data_miss,
        model_type = imputation_model,
        config = config,
        seed = rep_seeds$seed_brms
      )

      # 4. Extract imputed datasets
      imputed_datasets <- extract_imputed_datasets(brms_fit, data_miss, M = config$M)

      # 5. Fit lavaan.mi model to imputed data
      fit_mi <- fit_mediation_mi(imputed_datasets, model_type = "full_mediation")

      # 6. Compute RIV
      riv_result <- compute_RIV_from_mi(fit_mi)

      # 7. Compute all metrics
      metrics <- compute_study1_metrics(
        data_complete = data_complete,
        imputed_datasets = imputed_datasets,
        fit_mi = fit_mi,
        riv_result = riv_result,
        pop_params = pop_params,
        model_type = "full_mediation"
      )

      # Add condition identifiers
      metrics$rep_id <- rep_id
      metrics$mechanism <- mechanism
      metrics$imputation_model <- imputation_model

      # Add convergence info
      convergence <- check_brms_convergence(brms_fit)
      metrics$brms_converged <- convergence$converged
      metrics$brms_max_rhat <- convergence$max_rhat
    }

    # Add missingness info (common to both approaches)
    metrics$actual_prop_missing <- ampute_result$actual_prop_missing
    metrics$n_miss_M <- ampute_result$n_miss_M
    metrics$n_miss_Y <- ampute_result$n_miss_Y

    metrics

  }, error = function(e) {
    if (verbose) {
      cat(sprintf("  ERROR in rep %d: %s\n", rep_id, conditionMessage(e)))
    }
    list(
      rep_id = rep_id,
      mechanism = mechanism,
      imputation_model = imputation_model,
      error = conditionMessage(e),
      failed = TRUE
    )
  })

  return(result)
}


#' Run All Replications for One Condition (Serial)
#'
#' Runs all replications for a single mechanism × imputation model condition
#'
#' @param mechanism Character. "MCAR", "MAR", or "MNAR"
#' @param imputation_model Character. "matched", "saturated", or "uncongenial"
#' @param config List. Configuration
#' @param seeds List. Pre-generated seeds
#' @param verbose Logical. Print progress?
#' @return List of results for each replication
run_condition_serial <- function(mechanism, imputation_model, config, seeds, verbose = TRUE) {

  if (verbose) {
    cat(sprintf("\n=== Condition: %s × %s ===\n", mechanism, imputation_model))
  }

  results <- vector("list", config$n_reps)

  for (rep in seq_len(config$n_reps)) {
    if (verbose && rep %% 10 == 0) {
      cat(sprintf("  Progress: %d/%d\n", rep, config$n_reps))
    }

    results[[rep]] <- run_single_replication(
      rep_id = rep,
      mechanism = mechanism,
      imputation_model = imputation_model,
      config = config,
      seeds = seeds,
      verbose = FALSE
    )
  }

  # Count successes/failures
  n_failed <- sum(sapply(results, function(r) isTRUE(r$failed)))
  if (verbose) {
    cat(sprintf("  Completed: %d/%d (failed: %d)\n",
                config$n_reps - n_failed, config$n_reps, n_failed))
  }

  return(results)
}


#' Create Replication Task
#'
#' Creates a task specification for parallel execution
#'
#' @param rep_id Integer
#' @param mechanism Character
#' @param imputation_model Character
#' @return List with task specification
create_task <- function(rep_id, mechanism, imputation_model) {
  list(
    rep_id = rep_id,
    mechanism = mechanism,
    imputation_model = imputation_model
  )
}


#' Run Single Task (for parallel execution)
#'
#' Wrapper for run_single_replication that takes a task specification
#'
#' @param task List with rep_id, mechanism, imputation_model
#' @param config List. Configuration
#' @param seeds List. Pre-generated seeds
#' @return Result from run_single_replication()
run_task <- function(task, config, seeds) {
  run_single_replication(
    rep_id = task$rep_id,
    mechanism = task$mechanism,
    imputation_model = task$imputation_model,
    config = config,
    seeds = seeds,
    verbose = FALSE
  )
}


#' Generate All Tasks
#'
#' Creates list of all tasks for the full simulation
#'
#' @param config List. Configuration
#' @return List of task specifications
generate_all_tasks <- function(config) {

  tasks <- list()
  task_id <- 1

  for (rep in seq_len(config$n_reps)) {
    for (mech in config$mechanisms) {
      for (model in config$imputation_models) {
        tasks[[task_id]] <- create_task(rep, mech, model)
        task_id <- task_id + 1
      }
    }
  }

  return(tasks)
}


#' Organize Results by Condition
#'
#' Takes flat list of results and organizes by mechanism and imputation model
#'
#' @param results_flat List of results from parallel execution
#' @param config List. Configuration
#' @return Nested list: results[[mechanism]][[imputation_model]] = list of rep results
organize_results <- function(results_flat, config) {

  organized <- list()

  for (mech in config$mechanisms) {
    organized[[mech]] <- list()
    for (model in config$imputation_models) {
      organized[[mech]][[model]] <- list()
    }
  }

  for (result in results_flat) {
    if (!is.null(result) && !is.null(result$mechanism)) {
      mech <- result$mechanism
      model <- result$imputation_model
      rep_id <- result$rep_id

      organized[[mech]][[model]][[rep_id]] <- result
    }
  }

  return(organized)
}
