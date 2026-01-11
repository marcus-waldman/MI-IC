# ============================================================================
# Study 1: Pedagogical Simulation - Main Script
# ============================================================================
# Validates the bias formula: Empirical Bias ≈ 0.5 × tr(RIV)
# Demonstrates that congeniality (not ignorability) is the critical assumption
#
# Design: 3×3 factorial
#   - Mechanisms: MCAR, MAR, MNAR
#   - Imputation models: matched, saturated, uncongenial
# ============================================================================

# Clear environment
rm(list = ls())

# Set working directory to script location (if running interactively)
if (interactive()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# ============================================================================
# 1. Load Configuration and Source Functions
# ============================================================================

cat("=== Study 1: Pedagogical Simulation ===\n\n")
cat("Loading configuration and functions...\n")

# Load config
source("config.R")

# Source all utility functions
source("R/00_setup_seeds.R")
source("R/01_generate_data.R")
source("R/02_ampute_data.R")
source("R/03_impute_brms.R")
source("R/04_extract_imputations.R")
source("R/05_fit_lavaan_mi.R")
source("R/06_compute_metrics.R")
source("R/07_run_replication.R")

# Load required packages
required_packages <- c(
  "parallel",
  "pbapply",
  "brms",
  "lavaan",
  "lavaan.mi",
  "mice",
  "posterior",
  "SeedMaker"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing package '%s'...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

cat("  Packages loaded.\n")

# ============================================================================
# 2. Setup Seeds
# ============================================================================

cat("\nSetting up seeds...\n")
seeds <- setup_seeds(config)
print_seed_summary(seeds, config)

# ============================================================================
# 3. Run Complete-Data Validation (Optional but Recommended)
# ============================================================================

run_validation <- TRUE  # Set to FALSE to skip

if (run_validation) {
  cat("\n=== Complete-Data Validation ===\n")
  cat("Checking parameter recovery and CI coverage with complete data...\n")

  # Run validation with subset of replications
  n_validation_reps <- min(100, config$n_reps)

  validation_results <- run_complete_data_validation(
    n_reps = n_validation_reps,
    config = config,
    seeds = seeds,
    verbose = TRUE
  )

  # Summarize
  pop_params <- get_population_params(config)
  validation_summary <- summarize_validation(validation_results, pop_params)
  print_validation_summary(validation_summary)

  # Save validation results
  saveRDS(validation_results, "results/validation_results.rds")
  saveRDS(validation_summary, "results/validation_summary.rds")

  cat("\nValidation complete. Results saved to results/\n")
}

# ============================================================================
# 4. Run Main Simulation
# ============================================================================

cat("\n=== Main Simulation ===\n")
cat(sprintf("Design: %d mechanisms × %d imputation models × %d replications\n",
            length(config$mechanisms), length(config$imputation_models), config$n_reps))
cat(sprintf("Total tasks: %d\n",
            length(config$mechanisms) * length(config$imputation_models) * config$n_reps))

# Determine parallelization
n_cores <- get_n_cores(config)
cat(sprintf("Using %d cores for parallel processing\n\n", n_cores))

# Generate all tasks
tasks <- generate_all_tasks(config)
n_tasks <- length(tasks)

# Setup cluster
cat("Setting up parallel cluster...\n")
cl <- parallel::makeCluster(n_cores, type = config$parallel_type)

# Export required objects to workers
parallel::clusterExport(cl, c(
  "config", "seeds",
  # Functions
  "run_task", "run_single_replication",
  "get_replication_seeds", "get_population_params",
  "generate_mediation_data_from_config", "generate_mediation_data",
  "derive_mediation_covariance", "derive_mediation_covariance_from_config",
  "ampute_mediation",
  "fit_brms_imputation_from_config", "fit_brms_imputation", "build_imputation_formula",
  "extract_imputed_datasets",
  "fit_mediation_mi", "get_mediation_syntax", "compute_RIV_from_mi",
  "fit_mediation_complete", "check_parameter_coverage", "extract_lavaan_params",
  "get_pooled_estimates",
  "compute_study1_metrics", "get_loglik_per_imputation", "compute_loglik_at_params",
  "check_brms_convergence"
))

# Load packages on workers
parallel::clusterEvalQ(cl, {
  library(brms)
  library(lavaan)
  library(lavaan.mi)
  library(mice)
  library(posterior)
})

cat("Cluster ready.\n\n")

# Run simulation with progress bar
cat("Running simulation...\n")
start_time <- Sys.time()

results_flat <- pbapply::pblapply(
  X = tasks,
  FUN = function(task) run_task(task, config, seeds),
  cl = cl
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

# Stop cluster
parallel::stopCluster(cl)

cat(sprintf("\nSimulation complete. Elapsed time: %.1f minutes\n", as.numeric(elapsed)))

# ============================================================================
# 5. Organize and Summarize Results
# ============================================================================

cat("\n=== Organizing Results ===\n")

# Organize by condition
results_organized <- organize_results(results_flat, config)

# Create summary table
results_table <- create_results_table(results_organized, config)

# Check bias formula
results_table <- check_bias_formula(results_table)

# Check parameter recovery (quality check for pooled estimates from Rubin's rules)
results_table <- check_parameter_recovery(results_table)

# Print full summary (bias formula + parameter recovery)
print_full_results_summary(results_table)

# ============================================================================
# 6. Save Results
# ============================================================================

cat("\n=== Saving Results ===\n")

# Create results directory if needed
dir.create("results", showWarnings = FALSE)

# Save all results
saveRDS(results_flat, "results/results_flat.rds")
saveRDS(results_organized, "results/results_organized.rds")
saveRDS(results_table, "results/results_table.rds")

# Save summary as CSV for easy viewing
write.csv(results_table, "results/results_summary.csv", row.names = FALSE)

cat("Results saved to results/\n")
cat("  - results_flat.rds: Raw results list\n")
cat("  - results_organized.rds: Results organized by condition\n")
cat("  - results_table.rds: Summary table\n")
cat("  - results_summary.csv: Summary in CSV format\n")

# ============================================================================
# 7. Quick Summary
# ============================================================================

cat("\n=== Quick Summary ===\n")
cat(sprintf("Total replications: %d × %d conditions = %d tasks\n",
            config$n_reps, nrow(results_table), n_tasks))
cat(sprintf("Elapsed time: %.1f minutes (%.2f min/task)\n",
            as.numeric(elapsed), as.numeric(elapsed) / n_tasks))

# Key finding
congenial_results <- results_table[results_table$imputation_model != "uncongenial", ]
uncongenial_results <- results_table[results_table$imputation_model == "uncongenial", ]

cat("\nKey Findings:\n")
cat(sprintf("  Congenial conditions (matched, saturated): mean |rel diff| = %.1f%%\n",
            mean(abs(congenial_results$rel_diff), na.rm = TRUE) * 100))
cat(sprintf("  Uncongenial condition: mean |rel diff| = %.1f%%\n",
            mean(abs(uncongenial_results$rel_diff), na.rm = TRUE) * 100))

cat("\nDone!\n")
