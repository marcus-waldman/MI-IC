# ============================================================================
# Term Decomposition Verification Simulation (Analytic Version)
# ============================================================================
# Verify Taylor expansion terms for imputation bias in saturated MVN
# Uses lavaan for MLE and analytic Q-function (no imputation)
# ============================================================================

# Load packages
library(lavaan)
library(MASS)
library(pbapply)
library(parallel)
library(SeedMaker)

# Source analytic MVN functions
source("simulations/000-premise-very-basic-Y/mvn_analytics.R")
source("simulations/utils/generate_data.R")  # For impose_missingness with mice::ampute

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

N <- 2000        # Sample size (large to minimize variance)
p <- 5           # Dimension (Q = 20 parameters)
gamma <- 0.6     # Missing rate (high to maximize RIV)
n_reps <- 10000    # Number of replications
seed_base <- 12345

# True parameters
mu_star <- rep(0, p)

# Generate Sigma_star from Wishart distribution
set.seed(999)  # Fixed seed for reproducibility of true parameters
wishart_df <- p + 10
Sigma_star <- rWishart(1, df = wishart_df, Sigma = diag(p))[,,1] / wishart_df

# Variable names
var_names <- paste0("Y", 1:p)

# Generate saturated model syntax once
model_syntax <- generate_saturated_model(var_names)

# ============================================================================
# SINGLE REPLICATION FUNCTION
# ============================================================================

run_one_replication <- function(rep_idx, seed) {

  set.seed(seed)

  # -------------------------------------------------------------------------
  # Step 1: Generate complete data
  # -------------------------------------------------------------------------
  Y_complete <- mvrnorm(N, mu_star, Sigma_star)
  colnames(Y_complete) <- var_names

  # -------------------------------------------------------------------------
  # Step 2: Impose MCAR missingness (monotone pattern via mice::ampute)
  # -------------------------------------------------------------------------
  miss_result <- impose_missingness(
    data = Y_complete,
    missing_rate = gamma,
    pattern = "monotone",
    prop_complete = 0.4,
    seed = seed + 1000
  )
  Y_miss <- miss_result$data_miss
  colnames(Y_miss) <- var_names

  # -------------------------------------------------------------------------
  # Step 3: Fit lavaan models
  # -------------------------------------------------------------------------

  # Fit to complete data
  fit_com <- tryCatch({
    sem(model_syntax, data = as.data.frame(Y_complete),
        meanstructure = TRUE, likelihood = "normal")
  }, error = function(e) NULL)

  # Fit to incomplete data (FIML)
  fit_obs <- tryCatch({
    sem(model_syntax, data = as.data.frame(Y_miss),
        meanstructure = TRUE, missing = "fiml")
  }, error = function(e) NULL)

  if (is.null(fit_com) || is.null(fit_obs)) {
    warning(paste("Lavaan did not converge for replication", rep_idx))
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # Step 4: Extract MLE and Information Matrices
  # -------------------------------------------------------------------------

  # Get MLE from observed data fit (θ̂_obs)
  mle_obs <- get_mle_from_lavaan(fit_obs)
  mu_hat_obs <- mle_obs$mu
  Sigma_hat_obs <- mle_obs$Sigma
  theta_hat_obs <- mle_obs$theta_vec

  # Get MLE from complete data fit (θ̂_com)
  mle_com <- get_mle_from_lavaan(fit_com)
  mu_hat_com <- mle_com$mu
  Sigma_hat_com <- mle_com$Sigma
  theta_hat_com <- mle_com$theta_vec

  # True parameter vector
  theta_star <- c(mu_star, vech(Sigma_star))

  # Information matrices
  I_com <- get_info_matrix(fit_com)  # Complete data info at MLE
  I_obs <- get_info_matrix(fit_obs)  # Observed data info at MLE

  # -------------------------------------------------------------------------
  # Step 5: Compute Q-function and log-likelihood
  # -------------------------------------------------------------------------

  # Analytic Q-function at θ̂_obs
  Q_at_hat <- compute_Q_analytic(Y_miss, mu_hat_obs, Sigma_hat_obs)

  # Complete-data log-likelihood at θ̂_com
  ell_com_at_hat <- mvn_loglik(Y_complete, mu_hat_com, Sigma_hat_com)

  # Total Bias = Q(θ̂_obs) - ℓ_com(θ̂_com)
  total_bias <- Q_at_hat - ell_com_at_hat

  # -------------------------------------------------------------------------
  # Step 6: Compute Term 1 (at true parameters)
  # -------------------------------------------------------------------------
  Q_at_star <- compute_Q_analytic(Y_miss, mu_star, Sigma_star)
  ell_com_at_star <- mvn_loglik(Y_complete, mu_star, Sigma_star)
  term1 <- Q_at_star - ell_com_at_star

  # -------------------------------------------------------------------------
  # Step 7: Compute Term 2 (score difference)
  # -------------------------------------------------------------------------
  # Term 2 = (theta_hat - theta_star)' * (S_Q - S_com) evaluated at theta_star
  S_com <- mvn_score(Y_complete, mu_star, Sigma_star)

  # Compute S_Q using numerical differentiation (analytic version has issues)
  Q_len <- length(theta_star)
  S_Q <- numeric(Q_len)
  eps <- 1e-6

  for (j in 1:Q_len) {
    theta_plus <- theta_star
    theta_minus <- theta_star
    theta_plus[j] <- theta_plus[j] + eps
    theta_minus[j] <- theta_minus[j] - eps

    if (j <= p) {
      mu_plus <- theta_plus[1:p]
      mu_minus <- theta_minus[1:p]
      Sigma_plus <- Sigma_star
      Sigma_minus <- Sigma_star
    } else {
      mu_plus <- mu_star
      mu_minus <- mu_star
      Sigma_plus <- vech_inv(theta_plus[(p+1):Q_len], p)
      Sigma_minus <- vech_inv(theta_minus[(p+1):Q_len], p)
    }

    Q_plus <- compute_Q_analytic(Y_miss, mu_plus, Sigma_plus)
    Q_minus <- compute_Q_analytic(Y_miss, mu_minus, Sigma_minus)
    S_Q[j] <- (Q_plus - Q_minus) / (2 * eps)
  }

  term2 <- sum((theta_hat_obs - theta_star) * (S_Q - S_com))

  # -------------------------------------------------------------------------
  # Step 8: Compute Term 3 (Hessian difference)
  # -------------------------------------------------------------------------
  # Term 3 = -0.5 * (theta_hat_obs - theta_star)' * (J_Q - J_com) * (theta_hat_obs - theta_star)
  # For MVN, both J_Q and J_com equal the Fisher information at Sigma_star
  # (the Hessian doesn't depend on Y, only on Sigma)
  # So Term 3 = 0 analytically
  J_com <- mvn_fisher_info(N, Sigma_star)
  J_Q <- J_com  # Same for MVN since Hessian doesn't depend on data
  J_diff <- J_Q - J_com
  term3 <- -0.5 * as.numeric(t(theta_hat_obs - theta_star) %*% J_diff %*% (theta_hat_obs - theta_star))

  # -------------------------------------------------------------------------
  # Step 9: Taylor sum and residual
  # -------------------------------------------------------------------------
  taylor_sum <- term1 + term2 + term3
  residual <- total_bias - taylor_sum

  # -------------------------------------------------------------------------
  # Step 10: Compute tr(RIV) analytically
  # -------------------------------------------------------------------------
  tr_RIV <- compute_tr_RIV_analytic(I_com, I_obs)

  # -------------------------------------------------------------------------
  # Return results
  # -------------------------------------------------------------------------
  return(data.frame(
    replication = rep_idx,
    N = N,
    p = p,
    gamma = gamma,
    total_bias = total_bias,
    term1 = term1,
    term2 = term2,
    term3 = term3,
    taylor_sum = taylor_sum,
    residual = residual,
    Q_at_hat = Q_at_hat,
    Q_at_star = Q_at_star,
    ell_com_at_hat = ell_com_at_hat,
    ell_com_at_star = ell_com_at_star,
    tr_RIV = tr_RIV,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# RUN SIMULATION
# ============================================================================

cat("=== Term Decomposition Verification (Analytic) ===\n")
cat(sprintf("N=%d, p=%d, gamma=%.2f, n_reps=%d\n", N, p, gamma, n_reps))
cat(sprintf("Q = %d parameters (mu=%d, vech(Sigma)=%d)\n\n",
            p + p*(p+1)/2, p, p*(p+1)/2))

# Generate seeds
rep_seeds <- SeedMaker::seed_maker(
    seed = seed_base,
    level_names = "sim",
    n_per_level = n_reps
  )

if (n_reps == 1) {
  # Single replication - no parallelization needed
  cat("Running single replication...\n")
  t1 = proc.time()
  results_list <- list(run_one_replication(1, rep_seeds[[1]]))
  cat("\nRun time:\n")
  print(proc.time()-t1)
} else {
  # Multiple replications - use parallel processing
  cat("Setting up parallel processing...\n")
  n_cores <- min(25, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)

  # Export required objects and functions
  clusterExport(cl, c(
    "N", "p", "gamma", "mu_star", "Sigma_star", "var_names", "model_syntax",
    "vech", "vech_inv", "duplication_matrix",
    "mvn_loglik", "mvn_score", "mvn_fisher_info",
    "generate_saturated_model", "get_lavaan_to_theta_perm",
    "get_info_matrix", "get_mle_from_lavaan",
    "compute_Q_analytic", "compute_tr_RIV_analytic",
    "run_one_replication", "rep_seeds",
    "impose_missingness"
  ))

  clusterEvalQ(cl, {
    library(lavaan)
    library(MASS)
    library(mice)  # For impose_missingness -> ampute
  })

  cat(sprintf("Running %d replications on %d cores...\n", n_reps, n_cores))

  # Run with progress bar
  results_list <- pblapply(1:n_reps, function(i) {
    run_one_replication(i, rep_seeds[[i]])
  }, cl = cl)

  stopCluster(cl)
}

# Filter out failed replications
results_list <- results_list[!sapply(results_list, is.null)]
n_success <- length(results_list)

if (n_success == 0) {
  stop("All replications failed!")
}

cat(sprintf("\nSuccessful replications: %d / %d\n", n_success, n_reps))

# Combine results
results_df <- do.call(rbind, results_list)

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================

cat("\n=== Results ===\n\n")

if (n_reps == 1) {
  # Single replication - show raw values
  r <- results_df[1, ]

  cat("| Metric         | Value        | Ratio to tr(RIV) |\n")
  cat("|----------------|--------------|------------------|\n")
  cat(sprintf("| Total Bias     | %12.4f | %16.4f |\n", r$total_bias, r$total_bias / r$tr_RIV))
  cat(sprintf("| Term 1         | %12.4f | %16.4f |\n", r$term1, r$term1 / r$tr_RIV))
  cat(sprintf("| Term 2         | %12.4f | %16.4f |\n", r$term2, r$term2 / r$tr_RIV))
  cat(sprintf("| Term 3         | %12.4f | %16.4f |\n", r$term3, r$term3 / r$tr_RIV))
  cat(sprintf("| Taylor Sum     | %12.4f | %16.4f |\n", r$taylor_sum, r$taylor_sum / r$tr_RIV))
  cat(sprintf("| Residual       | %12.4f | %16.4f |\n", r$residual, r$residual / r$tr_RIV))
  cat(sprintf("| tr(RIV)        | %12.4f | %16.4f |\n", r$tr_RIV, 1.0))

} else {
  # Multiple replications - show means and SDs
  mean_tr_RIV <- mean(results_df$tr_RIV)

  summary_df <- data.frame(
    Metric = c("Total Bias", "Term 1", "Term 2", "Term 3", "Taylor Sum", "Residual", "tr(RIV)"),
    Mean = c(
      mean(results_df$total_bias),
      mean(results_df$term1),
      mean(results_df$term2),
      mean(results_df$term3),
      mean(results_df$taylor_sum),
      mean(results_df$residual),
      mean_tr_RIV
    ),
    SD = c(
      sd(results_df$total_bias),
      sd(results_df$term1),
      sd(results_df$term2),
      sd(results_df$term3),
      sd(results_df$taylor_sum),
      sd(results_df$residual),
      sd(results_df$tr_RIV)
    )
  )
  summary_df$Ratio_to_trRIV <- summary_df$Mean / mean_tr_RIV

  print(summary_df, digits = 4)
}

# Save results
dir.create("simulations/000-premise-very-basic-Y/results", showWarnings = FALSE)
saveRDS(results_df, "simulations/000-premise-very-basic-Y/results/term_decomposition_analytic.rds")
cat("\nResults saved to: simulations/000-premise-very-basic-Y/results/term_decomposition_analytic.rds\n")
