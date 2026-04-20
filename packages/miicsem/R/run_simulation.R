#' Generate Seeds for All Replications (3 stages)
#'
#' Produces three seeds per replication — one each for data generation,
#' amputation, and imputation — using \code{SeedMaker}.
#'
#' @param master_seed Master seed (integer).
#' @param n_reps Number of replications (integer).
#' @return List of \code{n_reps} elements, each with \code{data},
#'   \code{ampute}, \code{impute} integer seeds.
#' @export
generate_seeds <- function(master_seed, n_reps) {
  if (!requireNamespace("SeedMaker", quietly = TRUE)) {
    stop("Package 'SeedMaker' required. Install with: install.packages('SeedMaker')")
  }

  raw_seeds <- SeedMaker::seed_maker(
    seed        = master_seed,
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


#' Run the Full Study 2 SEM Simulation
#'
#' Primary user-facing entry point. Runs the Bollen et al. (2014) SIM1
#' simulation across a grid of sample sizes and missing rates, computing
#' all seven information criteria (AIC_com, BIC_com, AIC_adhoc, BIC_adhoc,
#' AICcd, MI-AIC, MI-BIC) for each of 12 candidate models. Uses
#' \code{parallel} + \code{pbapply} for multi-core execution with
#' per-condition checkpointing to \code{.rds} files.
#'
#' @param n_reps Number of replications per condition. Default 1000.
#' @param n_cores Number of parallel workers. Default:
#'   \code{parallel::detectCores()}.
#' @param seed Master seed for reproducibility. Default 32897891.
#' @param sample_sizes Numeric vector of sample sizes to evaluate.
#'   Default \code{c(100, 250, 500, 1000, 5000)}.
#' @param miss_rates Numeric vector of missing rates. Default
#'   \code{c(0.10, 0.25, 0.40)}.
#' @param M Number of imputations per replication. Default 20.
#' @param results_dir Directory for per-condition \code{.rds} checkpoints.
#'   If \code{NULL}, results are returned but not saved.
#' @param verbose Print progress messages? Default \code{TRUE}.
#'
#' @return Named list of condition results. Each element corresponds to
#'   one \code{n x miss_rate} combination and contains \code{n_reps}
#'   per-replication results (see \code{\link{run_one_rep}}).
#'
#' @examples
#' \dontrun{
#' # Quick pilot
#' res <- run_simulation(n_reps = 10,
#'                       sample_sizes = 250,
#'                       miss_rates = 0.25)
#'
#' # Full grid on HPC
#' res <- run_simulation(n_reps = 1000,
#'                       results_dir = "results/")
#' }
#' @export
run_simulation <- function(n_reps      = 1000L,
                           n_cores     = parallel::detectCores(),
                           seed        = 32897891L,
                           sample_sizes = c(100, 250, 500, 1000, 5000),
                           miss_rates   = c(0.10, 0.25, 0.40),
                           M           = 20L,
                           results_dir = NULL,
                           verbose     = TRUE) {

  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Package 'pbapply' required.")
  }

  config <- get_config(
    n_reps       = n_reps,
    master_seed  = seed,
    sample_sizes = sample_sizes,
    miss_rates   = miss_rates,
    M            = M
  )

  if (!is.null(results_dir) && !dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  if (verbose) message("Generating seeds...")
  seeds <- generate_seeds(config$master_seed, config$n_reps)

  conditions <- expand.grid(
    n         = config$sample_sizes,
    miss_rate = config$miss_rates,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message(sprintf("Conditions: %d (%d sample sizes x %d miss rates)",
                    nrow(conditions), length(config$sample_sizes),
                    length(config$miss_rates)))
    message(sprintf("Using %d cores", n_cores))
  }

  all_results <- list()

  for (cond_idx in seq_len(nrow(conditions))) {
    n_val  <- conditions$n[cond_idx]
    mr_val <- conditions$miss_rate[cond_idx]
    cond_label <- sprintf("n=%d_mr=%.2f", n_val, mr_val)

    if (verbose) {
      message(sprintf("\n=== Condition %d/%d: %s ===",
                      cond_idx, nrow(conditions), cond_label))
    }

    if (!is.null(results_dir)) {
      rds_file <- file.path(results_dir, sprintf("results_%s.rds", cond_label))
      if (file.exists(rds_file)) {
        if (verbose) message("  Loading cached results...")
        all_results[[cond_label]] <- readRDS(rds_file)
        next
      }
    }

    cl <- parallel::makeCluster(n_cores)

    parallel::clusterExport(cl,
      varlist = c("seeds", "config", "n_val", "mr_val"),
      envir   = environment())

    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        requireNamespace("miicsem", quietly = TRUE)
        library(lavaan)
        library(mice)
        library(MASS)
      })
    })

    cond_results <- tryCatch({
      pbapply::pblapply(seq_len(config$n_reps), function(rep) {
        miicsem::run_one_rep(
          rep_id      = rep,
          n           = n_val,
          miss_rate   = mr_val,
          config      = config,
          seed_data   = seeds[[rep]]$data,
          seed_ampute = seeds[[rep]]$ampute,
          seed_impute = seeds[[rep]]$impute
        )
      }, cl = cl)
    }, finally = {
      parallel::stopCluster(cl)
    })

    if (!is.null(results_dir)) {
      saveRDS(cond_results, rds_file)
      if (verbose) message(sprintf("  Saved to %s", rds_file))
    }

    all_results[[cond_label]] <- cond_results
  }

  if (!is.null(results_dir)) {
    combined_file <- file.path(results_dir, "results_combined.rds")
    saveRDS(all_results, combined_file)
    if (verbose) message(sprintf("\nAll results saved to %s", combined_file))
  }

  return(all_results)
}
