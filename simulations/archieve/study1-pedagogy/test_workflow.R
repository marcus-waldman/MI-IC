# ============================================================================
# Test Script: Validate Study 1 Workflow
# ============================================================================
# Runs a minimal test (1 replication, 1 condition) to verify all components
# work together before running the full simulation.
# ============================================================================

# Clear environment
rm(list = ls())

# Set working directory to script location (if running interactively)
if (interactive()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

cat("=== Study 1: Workflow Validation Test ===\n\n")

# ============================================================================
# 1. Load Functions
# ============================================================================

cat("1. Loading functions...\n")

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
required_packages <- c("brms", "lavaan", "lavaan.mi", "mice", "posterior", "SeedMaker")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("   Installing package '%s'...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("   All packages loaded.\n\n")

# ============================================================================
# 2. Create Test Configuration (minimal settings)
# ============================================================================

cat("2. Creating test configuration...\n")

test_config <- list(
  # Small sample for testing
  n = 100,
  prop_missing = 0.3,
  n_reps = 1,  # Just 1 replication for testing

  # Use 1 core for testing
  n_cores = 1,
  parallel_type = "PSOCK",

  # Seeds
  master_seed = 12345,
  seeds_file = "results/test_seeds.rds",

  # Full brms settings for convergence check
  brms_backend = "rstan",
  brms_chains = 1,
  brms_warmup = 1000,   # Proper warmup
  brms_iter = 2000,     # 1000 post-warmup samples
  brms_refresh = 0,

  # Population parameters
  beta_xm = 0.5,
  beta_my = 0.5,
  beta_xy = 0,
  sigma2_m = 1,
  sigma2_y = 1,
  var_x = 1,

  # Test single condition
  mechanisms = c("MCAR"),
  imputation_models = c("matched"),

  # Seed offsets
  mechanism_offsets = c(MCAR = 0, MAR = 1000, MNAR = 2000),
  model_offsets = c(matched = 0, saturated = 100, uncongenial = 200)
)

cat("   Configuration:\n")
cat(sprintf("     n = %d, prop_missing = %.1f%%\n", test_config$n, test_config$prop_missing * 100))
cat(sprintf("     brms: %d warmup, %d iter (= %d samples)\n",
            test_config$brms_warmup, test_config$brms_iter,
            test_config$brms_iter - test_config$brms_warmup))
cat("     This will take ~3-5 minutes for proper convergence check\n")
cat("\n")

# ============================================================================
# 3. Test Data Generation
# ============================================================================

cat("3. Testing data generation...\n")

set.seed(test_config$master_seed)
data_complete <- generate_mediation_data_from_config(test_config, seed = 101)

cat(sprintf("   Generated data: n = %d\n", nrow(data_complete)))
cat(sprintf("   Column names: %s\n", paste(names(data_complete), collapse = ", ")))
cat(sprintf("   Sample correlations: cor(X,M) = %.3f, cor(M,Y) = %.3f\n",
            cor(data_complete$X, data_complete$M),
            cor(data_complete$M, data_complete$Y)))
cat("\n")

# ============================================================================
# 4. Test Amputation
# ============================================================================

cat("4. Testing amputation (MCAR)...\n")

ampute_result <- ampute_mediation(
  data = data_complete,
  prop_missing = test_config$prop_missing,
  mechanism = "MCAR",
  seed = 201
)

data_miss <- ampute_result$data_miss

cat(sprintf("   Missing M: %d (%.1f%%)\n",
            ampute_result$n_miss_M,
            100 * ampute_result$n_miss_M / test_config$n))
cat(sprintf("   Missing Y: %d (%.1f%%)\n",
            ampute_result$n_miss_Y,
            100 * ampute_result$n_miss_Y / test_config$n))
cat("\n")

# ============================================================================
# 5. Test brms Imputation
# ============================================================================

cat("5. Testing brms imputation (matched model)...\n")
cat("   This will take ~3-5 minutes with 1000 warmup + 1000 sampling iterations...\n")

brms_fit <- fit_brms_imputation_from_config(
  data_miss = data_miss,
  model_type = "matched",
  config = test_config,
  seed = 301
)

# Check convergence
convergence <- check_brms_convergence(brms_fit)
cat(sprintf("   Converged: %s (max Rhat = %.3f)\n",
            ifelse(convergence$converged, "YES", "NO"),
            convergence$max_rhat))
cat("\n")

# ============================================================================
# 6. Test Imputation Extraction
# ============================================================================

cat("6. Testing imputation extraction...\n")

imputed_datasets <- extract_imputed_datasets(brms_fit, data_miss, M = NULL)

cat(sprintf("   Extracted %d imputed datasets\n", length(imputed_datasets)))

# Check first dataset has no missing
n_missing_first <- sum(is.na(imputed_datasets[[1]]))
cat(sprintf("   Missing values in first imputed dataset: %d\n", n_missing_first))
cat("\n")

# ============================================================================
# 7. Test lavaan.mi Analysis
# ============================================================================

cat("7. Testing lavaan.mi analysis...\n")

fit_mi <- fit_mediation_mi(imputed_datasets, model_type = "full_mediation")

cat("   Model fitted successfully.\n")

# Compute RIV
riv_result <- compute_RIV_from_mi(fit_mi)
cat(sprintf("   tr(RIV) = %.4f\n", riv_result$tr_RIV))

# Get pooled parameter estimates (Rubin's rules)
pop_params <- get_population_params(test_config)
pooled_params <- get_pooled_estimates(fit_mi)
pooled_params <- check_parameter_coverage(pooled_params, pop_params)

cat("\n   Pooled Parameter Estimates (Rubin's Rules via sem.mi):\n")
for (param_label in c("a", "b", "indirect")) {
  param_row <- pooled_params[pooled_params$label == param_label, ]
  if (nrow(param_row) > 0) {
    true_val <- param_row$true_value[1]
    est_val <- param_row$est[1]
    se_val <- param_row$se[1]
    bias <- param_row$bias[1]
    covered <- param_row$covered[1]
    cat(sprintf("     %s: est=%.3f (true=%.3f), SE=%.3f, bias=%+.3f, 95%% CI covers: %s\n",
                param_label, est_val, true_val, se_val, bias,
                ifelse(covered, "YES", "NO")))
  }
}
cat("\n   NOTE: For MCAR + matched model, estimates should be unbiased.\n")
cat("\n")

# ============================================================================
# 8. Test Complete-Data Analysis (for comparison)
# ============================================================================

cat("8. Testing complete-data analysis...\n")

fit_complete <- fit_mediation_complete(data_complete, model_type = "full_mediation")
ell_com <- lavaan::fitMeasures(fit_complete, "logl")

cat(sprintf("   Complete-data log-likelihood: %.4f\n", ell_com))
cat("\n")

# ============================================================================
# 9. Test Metrics Computation
# ============================================================================

cat("9. Testing metrics computation...\n")

# pop_params already defined in step 7

metrics <- compute_study1_metrics(
  data_complete = data_complete,
  imputed_datasets = imputed_datasets,
  fit_mi = fit_mi,
  riv_result = riv_result,
  pop_params = pop_params,
  model_type = "full_mediation"
)

cat(sprintf("   Q_bar (avg imputed loglik): %.4f\n", metrics$Q_bar))
cat(sprintf("   ell_com (complete loglik): %.4f\n", metrics$ell_com))
cat(sprintf("   Empirical bias (Q_bar - ell_com): %.4f\n", metrics$empirical_bias))
cat(sprintf("   Theoretical bias (0.5 * tr(RIV)): %.4f\n", metrics$theoretical_bias))
cat(sprintf("   Difference: %.4f\n", metrics$bias_difference))
cat(sprintf("   Ratio (empirical/theoretical): %.3f\n", metrics$bias_ratio))
cat("\n")

# ============================================================================
# 10. Test Diagnostic Decomposition (Levels 1-3)
# ============================================================================

cat("10. Testing bias decomposition diagnostics...\n\n")

cat("=== Level 1: Bias Decomposition (Approach A) ===\n")
cat(sprintf("   ell_com_at_pooled:         %.4f\n", metrics$ell_com_at_pooled))
cat(sprintf("   term_A1 (Imputation Bias): %.4f\n", metrics$term_A1))
cat(sprintf("   term_A2 (Estimation Mism): %.4f\n", metrics$term_A2))
cat(sprintf("   Sum (= empirical bias):    %.4f %s\n",
            metrics$term_A1 + metrics$term_A2,
            ifelse(abs((metrics$term_A1 + metrics$term_A2) - metrics$empirical_bias) < 1e-6, "✓", "✗")))
cat("\n")
cat("   Ratios (should ≈ 1.0 in expectation):\n")
cat(sprintf("   term_A1 / tr(RIV):         %.3f\n", metrics$term_A1_ratio))
cat(sprintf("   term_A2 / (-0.5×tr(RIV)):  %.3f\n", metrics$term_A2_ratio))
cat("\n")

cat("=== Level 2: Taylor Approximation Check ===\n")
cat(sprintf("   term_A2 (empirical):       %.4f\n", metrics$term_A2))
cat(sprintf("   term_A2_taylor (quadratic):%.4f\n", metrics$term_A2_taylor))
if (!is.na(metrics$term_A2_taylor)) {
  cat(sprintf("   Difference:                %.4f\n", metrics$term_A2 - metrics$term_A2_taylor))
}
cat("\n")

cat("=== Level 3: Information Matrices and Bilinear Form ===\n")
cat(sprintf("   tr(I_com):                 %.4f\n", metrics$tr_I_com))
cat(sprintf("   tr(I_obs):                 %.4f\n", metrics$tr_I_obs))
cat(sprintf("   tr(I_mis|obs):             %.4f\n", metrics$tr_I_mis_obs))
cat(sprintf("   term1_bilinear:            %.4f\n", metrics$term1_bilinear))
if (!is.na(metrics$term1_bilinear)) {
  cat(sprintf("   Compare to term_A1:        %.4f (diff: %.4f)\n",
              metrics$term_A1, metrics$term1_bilinear - metrics$term_A1))
}
cat("\n")

# ============================================================================
# 11. Summary
# ============================================================================

cat("=== Validation Test Summary ===\n\n")

checks <- c(
  "Data generation" = nrow(data_complete) == test_config$n,
  "Amputation" = ampute_result$n_miss_M > 0 && ampute_result$n_miss_Y > 0,
  "brms fit" = !is.null(brms_fit),
  "Imputation extraction" = length(imputed_datasets) > 0 && n_missing_first == 0,
  "lavaan.mi fit" = !is.null(fit_mi),
  "Pooled estimates (Rubin)" = nrow(pooled_params) > 0 && !is.na(pooled_params$est[1]),
  "RIV computation" = !is.na(riv_result$tr_RIV),
  "Metrics computation" = !is.na(metrics$empirical_bias),
  "Level 1 diagnostics" = !is.na(metrics$term_A1) && !is.na(metrics$term_A2),
  "Level 2 diagnostics" = !is.na(metrics$term_A2_taylor),
  "Level 3 diagnostics" = !is.na(metrics$term1_bilinear)
)

for (check_name in names(checks)) {
  status <- ifelse(checks[[check_name]], "PASS", "FAIL")
  cat(sprintf("  [%s] %s\n", status, check_name))
}

all_passed <- all(checks)
cat(sprintf("\nOverall: %s\n", ifelse(all_passed, "ALL TESTS PASSED", "SOME TESTS FAILED")))

if (all_passed) {
  cat("\nThe workflow is ready for full simulation.\n")
  cat("Run 'source(\"run_study1.R\")' to execute the full simulation.\n")
}
