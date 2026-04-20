# ============================================================================
# Parallel Simulation Driver for Study 2 (SEM)
# ============================================================================

# Source all modules
source_all <- function(sim_dir = NULL) {
  if (is.null(sim_dir)) {
    sim_dir <- "simulations/study2-bollen-sem/R"
  }
  source(file.path(sim_dir, "00_config.R"))
  source(file.path(sim_dir, "01_model_specifications.R"))
  source(file.path(sim_dir, "02_generate_and_ampute.R"))
  source(file.path(sim_dir, "03_fit_and_pool.R"))
  source(file.path(sim_dir, "04_compute_ic.R"))
  source(file.path(sim_dir, "05_run_replication.R"))
}


#' Generate Seeds for All Replications (3 stages)
#'
#' @param master_seed Master seed
#' @param n_reps Number of replications
#' @return List of n_reps elements, each with $data, $ampute, $impute seeds
generate_seeds <- function(master_seed, n_reps) {
  if (!requireNamespace("SeedMaker", quietly = TRUE)) {
    stop("Package 'SeedMaker' required. Install with: install.packages('SeedMaker')")
  }

  raw_seeds <- SeedMaker::seed_maker(
    seed = master_seed,
    level_names = c("rep", "stage"),
    n_per_level = c(n_reps, 3)
  )

  seeds <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    rep_name <- paste0("rep", i)
    seeds[[i]] <- list(
      data   = raw_seeds[[rep_name]][["stage1"]],
      ampute = raw_seeds[[rep_name]][["stage2"]],
      impute = raw_seeds[[rep_name]][["stage3"]]
    )
  }
  return(seeds)
}


#' Run Full Simulation
#'
#' @param config Configuration list from get_config()
#' @param n_cores Number of cores (NULL = auto)
#' @param results_dir Directory for saving intermediate results
run_simulation <- function(config, n_cores = NULL, results_dir = NULL) {
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Package 'pbapply' required.")
  }

  if (is.null(results_dir)) {
    results_dir <- "simulations/study2-bollen-sem/results"
  }
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # Generate seeds
  cat("Generating seeds...\n")
  seeds <- generate_seeds(config$master_seed, config$n_reps)

  # Condition grid
  conditions <- expand.grid(
    n         = config$sample_sizes,
    miss_rate = config$miss_rates,
    stringsAsFactors = FALSE
  )
  cat(sprintf("Conditions: %d (%d sample sizes x %d miss rates)\n",
              nrow(conditions), length(config$sample_sizes),
              length(config$miss_rates)))

  # Set up cluster
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  cat(sprintf("Using %d cores\n", n_cores))

  # Functions to export to workers
  export_fns <- c(
    "run_one_rep", "get_sim1_models", "get_config",
    "generate_complete_data", "ampute_data",
    "fit_complete", "fit_single_model_mi", "pool_mi",
    "eval_loglik_at_pooled", "fit_mi_models",
    "compute_all_ic", "compute_all_models_ic", "select_models"
  )

  # Run each condition
  all_results <- list()

  for (cond_idx in seq_len(nrow(conditions))) {
    n_val <- conditions$n[cond_idx]
    mr_val <- conditions$miss_rate[cond_idx]
    cond_label <- sprintf("n=%d_mr=%.2f", n_val, mr_val)

    cat(sprintf("\n=== Condition %d/%d: %s ===\n",
                cond_idx, nrow(conditions), cond_label))

    # Check for existing results (crash resilience)
    rds_file <- file.path(results_dir, sprintf("results_%s.rds", cond_label))
    if (file.exists(rds_file)) {
      cat("  Loading cached results...\n")
      all_results[[cond_label]] <- readRDS(rds_file)
      next
    }

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export to workers
    parallel::clusterExport(cl, export_fns, envir = environment())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(lavaan)
        library(mice)
        library(MASS)
      })
    })

    # Run replications with progress bar
    cond_results <- pbapply::pblapply(seq_len(config$n_reps), function(rep) {
      run_one_rep(
        rep_id      = rep,
        n           = n_val,
        miss_rate   = mr_val,
        config      = config,
        seed_data   = seeds[[rep]]$data,
        seed_ampute = seeds[[rep]]$ampute,
        seed_impute = seeds[[rep]]$impute
      )
    }, cl = cl)

    parallel::stopCluster(cl)
    on.exit(NULL) # Clear the on.exit since we already stopped

    # Save intermediate results
    saveRDS(cond_results, rds_file)
    cat(sprintf("  Saved to %s\n", rds_file))

    all_results[[cond_label]] <- cond_results
  }

  # Save combined results
  combined_file <- file.path(results_dir, "results_combined.rds")
  saveRDS(all_results, combined_file)
  cat(sprintf("\nAll results saved to %s\n", combined_file))

  return(all_results)
}


# --- Main execution ---
if (sys.nframe() == 0) {
  source_all()
  config <- get_config()
  run_simulation(config)
}
