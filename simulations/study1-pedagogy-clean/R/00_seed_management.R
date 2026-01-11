# ============================================================================
# Seed Management for Study 1 (Clean Implementation)
# ============================================================================
# Uses SeedMaker package to pre-generate reproducible seeds
# Simplified for analytical approach (no brms stage needed)
# ============================================================================

#' Generate Seeds Using SeedMaker
#'
#' @param master_seed Integer. Master seed for reproducibility
#' @param n_reps Integer. Number of replications
#' @return List with data and ampute seeds for each replication
generate_seeds <- function(master_seed, n_reps) {

  if (!requireNamespace("SeedMaker", quietly = TRUE)) {
    stop("Package 'SeedMaker' required. Install with: install.packages('SeedMaker')")
  }

  # Generate seeds with SeedMaker
  # 2 stages per replication: data generation, amputation
  raw_seeds <- SeedMaker::seed_maker(
    seed = master_seed,
    level_names = c("rep", "stage"),
    n_per_level = c(n_reps, 2)
  )

  # Reorganize into usable structure
  seeds <- vector("list", n_reps)
  names(seeds) <- sprintf("rep_%03d", 1:n_reps)

  for (i in 1:n_reps) {
    rep_name <- paste0("rep", i)
    seeds[[i]] <- list(
      data = raw_seeds[[rep_name]][["stage1"]],
      ampute = raw_seeds[[rep_name]][["stage2"]]
    )
  }

  return(seeds)
}


#' Print Seed Summary
#'
#' @param seeds List. Seeds from generate_seeds()
#' @param master_seed Integer. Master seed used
print_seed_summary <- function(seeds, master_seed) {
  cat("=== Seed Summary ===\n")
  cat(sprintf("Master seed: %d\n", master_seed))
  cat(sprintf("Replications: %d\n", length(seeds)))
  cat(sprintf("Stages per rep: 2 (data, ampute)\n"))
  cat("\nFirst 5 replications:\n")
  for (i in 1:min(5, length(seeds))) {
    cat(sprintf("  rep_%03d: data=%d, ampute=%d\n",
                i, seeds[[i]]$data, seeds[[i]]$ampute))
  }
}
