# ============================================================================
# Seed Management for Study 1
# ============================================================================
# Uses SeedMaker package to pre-generate reproducible seeds
# Seeds are saved to RDS for persistence across sessions
# ============================================================================

#' Setup Seeds for Simulation
#'
#' Loads seeds from file if exists, otherwise generates and saves them
#'
#' @param config List. Configuration from config.R
#' @param force_regenerate Logical. If TRUE, regenerate even if file exists
#' @return List of seeds with hierarchical structure
#' @details
#'   Seed structure (stage level + offsets):
#'   - Stage 1: Data generation (same across conditions within rep)
#'   - Stage 2: Amputation (use with mechanism offset)
#'   - Stage 3: brms imputation (use with model offset)
#'
#'   To reproduce specific replication:
#'   - seeds$rep_001$data, seeds$rep_001$ampute, seeds$rep_001$brms
#'
#' @examples
#'   source("config.R")
#'   seeds <- setup_seeds(config)
setup_seeds <- function(config, force_regenerate = FALSE) {

 if (!requireNamespace("SeedMaker", quietly = TRUE)) {
    stop("Package 'SeedMaker' required. Install with: install.packages('SeedMaker')")
  }

  seeds_path <- file.path(
    dirname(sys.frame(1)$ofile %||% "."),
    "..",
    config$seeds_file
  )

 # Normalize path
  seeds_path <- normalizePath(seeds_path, mustWork = FALSE)

  # Check if seeds file exists
 if (file.exists(seeds_path) && !force_regenerate) {
    message("Loading pre-generated seeds from: ", seeds_path)
    seeds <- readRDS(seeds_path)

    # Validate seed count matches config
    if (length(seeds) != config$n_reps) {
      warning(sprintf(
        "Seed file has %d replications but config specifies %d. Regenerating.",
        length(seeds), config$n_reps
      ))
      force_regenerate <- TRUE
    }
  }

  # Generate seeds if needed
  if (!file.exists(seeds_path) || force_regenerate) {
    message("Generating seeds with master seed: ", config$master_seed)

    seeds <- generate_seeds_seedmaker(
      master_seed = config$master_seed,
      n_reps = config$n_reps
    )

    # Ensure directory exists
    dir.create(dirname(seeds_path), recursive = TRUE, showWarnings = FALSE)

    # Save seeds
    saveRDS(seeds, seeds_path)
    message("Seeds saved to: ", seeds_path)
  }

  return(seeds)
}


#' Generate Seeds Using SeedMaker
#'
#' Internal function to generate hierarchical seeds
#'
#' @param master_seed Integer. Master seed for reproducibility
#' @param n_reps Integer. Number of replications
#' @return Named list of seeds
generate_seeds_seedmaker <- function(master_seed, n_reps) {

  # Generate seeds with SeedMaker
  # Level: replication, with 3 stages per replication
  raw_seeds <- SeedMaker::seed_maker(
    seed = master_seed,
    level_names = c("rep", "stage"),
    n_per_level = c(n_reps, 3)
  )

  # Reorganize into more usable structure
  seeds <- vector("list", n_reps)
  names(seeds) <- sprintf("rep_%03d", 1:n_reps)

  for (i in 1:n_reps) {
    rep_name <- paste0("rep", i)
    seeds[[i]] <- list(
      data = raw_seeds[[rep_name]][["stage1"]],
      ampute = raw_seeds[[rep_name]][["stage2"]],
      brms = raw_seeds[[rep_name]][["stage3"]]
    )
  }

  return(seeds)
}


#' Get Seeds for a Single Replication
#'
#' Extracts seeds for a specific replication with condition offsets applied
#'
#' @param seeds List. Full seeds object from setup_seeds()
#' @param rep_id Integer. Replication number (1 to n_reps)
#' @param mechanism Character. "MCAR", "MAR", or "MNAR"
#' @param imputation_model Character. "matched", "saturated", or "uncongenial"
#' @param config List. Configuration with offset definitions
#' @return List with seed_data, seed_ampute, seed_brms (with offsets applied)
#' @examples
#'   rep_seeds <- get_replication_seeds(seeds, rep_id = 1,
#'                                       mechanism = "MAR",
#'                                       imputation_model = "saturated",
#'                                       config = config)
get_replication_seeds <- function(seeds, rep_id, mechanism, imputation_model, config) {

  rep_name <- sprintf("rep_%03d", rep_id)

  if (!rep_name %in% names(seeds)) {
    stop(sprintf("Replication %d not found in seeds", rep_id))
  }

  base_seeds <- seeds[[rep_name]]

  # Apply offsets for condition-specific randomness
  mech_offset <- config$mechanism_offsets[[mechanism]]
  model_offset <- config$model_offsets[[imputation_model]]

  if (is.null(mech_offset)) {
    stop(sprintf("Unknown mechanism: %s", mechanism))
  }
  if (is.null(model_offset)) {
    stop(sprintf("Unknown imputation_model: %s", imputation_model))
  }

  return(list(
    seed_data = base_seeds$data,                              # Same across conditions
    seed_ampute = base_seeds$ampute + mech_offset,            # Varies by mechanism
    seed_brms = base_seeds$brms + mech_offset + model_offset  # Varies by full condition
  ))
}


#' Print Seed Summary
#'
#' Displays summary information about the seeds
#'
#' @param seeds List. Seeds object from setup_seeds()
#' @param config List. Configuration
print_seed_summary <- function(seeds, config) {
  cat("=== Seed Summary ===\n")
  cat(sprintf("Master seed: %d\n", config$master_seed))
  cat(sprintf("Replications: %d\n", length(seeds)))
  cat(sprintf("Stages per rep: %d (data, ampute, brms)\n", 3))
  cat("\nFirst 3 replications:\
")
  for (i in 1:min(3, length(seeds))) {
    cat(sprintf("  rep_%03d: data=%d, ampute=%d, brms=%d\n",
                i, seeds[[i]]$data, seeds[[i]]$ampute, seeds[[i]]$brms))
  }
  cat("\nMechanism offsets:", paste(names(config$mechanism_offsets),
                                    config$mechanism_offsets, sep = "=", collapse = ", "), "\n")
  cat("Model offsets:", paste(names(config$model_offsets),
                              config$model_offsets, sep = "=", collapse = ", "), "\n")
}
