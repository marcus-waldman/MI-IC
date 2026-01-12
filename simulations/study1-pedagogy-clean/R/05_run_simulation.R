# ============================================================================
# Run Study 1 Simulation (Parallel with pbapply)
# ============================================================================

#' Run One Replication with Separate Seeds
#'
#' @param seed_data Seed for data generation
#' @param seed_ampute Seed for amputation
#' @param n Sample size
#' @param mechanism "MCAR", "MAR", or "MNAR"
#' @param gamma Strength parameter for MAR/MNAR (default: 0.5)
#' @param cut1 Quantile for Y-only missing (default: 0.5)
#' @param cut2 Quantile for M+Y missing (default: 0.75)
#' @param true_params List with a, b, sigma2_M, sigma2_Y
#' @param imputation_model One of "A", "B", "C" (default: "A")
#' @return Named vector with results
run_one_rep <- function(seed_data, seed_ampute, n, mechanism, gamma, cut1, cut2, true_params,
                        imputation_model = "A") {

  # Generate data with data seed
  set.seed(seed_data)
  X <- rnorm(n)
  M <- true_params$a * X + rnorm(n, 0, sqrt(true_params$sigma2_M))
  Y <- true_params$b * M + rnorm(n, 0, sqrt(true_params$sigma2_Y))
  data_com <- data.frame(X = X, M = M, Y = Y)

  # Ampute data (monotone missingness only)
  data_mis <- ampute_data(data_com, mechanism = mechanism, gamma = gamma,
                          cut1 = cut1, cut2 = cut2, seed = seed_ampute)

  # Fit models
  theta_start <- c(0.5, 0.5, log(1), log(1))

  fit_c <- tryCatch(
    optim(theta_start, negloglik_complete, data = data_com, method = "BFGS", hessian = TRUE),
    error = function(e) NULL
  )
  fit_o <- tryCatch(
    optim(theta_start, negloglik_observed, data = data_mis, method = "BFGS", hessian = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit_c) || is.null(fit_o)) return(NULL)

  # Transform parameters (log-scale for variances)
  theta_c <- fit_c$par
  theta_c[3:4] <- exp(theta_c[3:4])

  theta_o <- fit_o$par
  theta_o[3:4] <- exp(theta_o[3:4])

  # Compute posterior covariance V_theta = I_obs^{-1}
  # Note: Hessian is of negative log-likelihood, so it's already the information matrix
  V_theta <- tryCatch(
    solve(fit_o$hessian),
    error = function(e) NULL
  )

  if (is.null(V_theta)) return(NULL)

  # Transform V_theta to original scale (delta method for exp transform)
  # For variance parameters, the Jacobian is diag(1, 1, exp(log_s2_M), exp(log_s2_Y))
  J <- diag(c(1, 1, theta_o[3], theta_o[4]))
  V_theta_orig <- J %*% V_theta %*% t(J)

  # Estimate imputation model parameters
  theta_impute <- estimate_imputation_params(data_mis, theta_o, model = imputation_model)

  # Estimate posterior covariance for imputation parameters
  V_theta_impute <- estimate_V_theta(data_mis, theta_impute, V_theta_orig, model = imputation_model)

  # Compute Q-function (both improper and proper)
  Q_results <- compute_Q_analytical_proper(data_mis, theta_o, theta_impute, V_theta_impute, model = imputation_model)

  # Complete-data log-likelihoods
  ell_com_at_com <- -fit_c$value
  ell_com_at_obs <- -negloglik_complete(c(theta_o[1], theta_o[2], log(theta_o[3]), log(theta_o[4])), data_com)

  # Bias decomposition (improper MI)
  term1_improper <- Q_results$Q_improper - ell_com_at_obs
  term2 <- ell_com_at_obs - ell_com_at_com
  total_improper <- term1_improper + term2

  # Bias decomposition (proper MI)
  term1_proper <- Q_results$Q_proper - ell_com_at_obs
  total_proper <- term1_proper + term2

  # Compute tr(RIV) from information matrices
  I_c <- fit_c$hessian
  I_o <- fit_o$hessian
  tr_RIV <- tryCatch(
    sum(diag(I_c %*% solve(I_o))) - 4,
    error = function(e) NA
  )

  c(
    # Improper MI results
    term1_improper = term1_improper,
    term2 = term2,
    total_improper = total_improper,

    # Proper MI results
    term1_proper = term1_proper,
    total_proper = total_proper,
    hessian_correction = Q_results$hessian_correction,

    # Reference quantities
    tr_RIV = tr_RIV,
    Q_improper = Q_results$Q_improper,
    Q_proper = Q_results$Q_proper,
    ell_com_at_obs = ell_com_at_obs,
    ell_com_at_com = ell_com_at_com,

    # Parameter estimates
    a_est = theta_o[1],
    b_est = theta_o[2],
    a_bias = theta_o[1] - true_params$a,
    b_bias = theta_o[2] - true_params$b
  )
}

#' Run Multiple Replications with Precomputed Seeds (Parallel)
#'
#' @param seeds List of seeds from generate_seeds() - each element has $data and $ampute
#' @param n Sample size per replication
#' @param mechanism "MCAR", "MAR", or "MNAR" (default: "MCAR")
#' @param gamma Strength parameter for MAR/MNAR (default: 0.5)
#' @param cut1 Quantile for Y-only missing (default: 0.5)
#' @param cut2 Quantile for M+Y missing (default: 0.75)
#' @param true_params List with a, b, sigma2_M, sigma2_Y
#' @param imputation_model One of "A", "B", "C" (default: "A")
#' @param n_cores Number of cores for parallel processing (default: auto-detect)
#' @return Data frame with results
run_simulation <- function(seeds, n, mechanism = "MCAR", gamma = 0.5,
                           cut1 = 0.5, cut2 = 0.75, true_params,
                           imputation_model = "A", n_cores = NULL) {

  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Package 'pbapply' required. Install with: install.packages('pbapply')")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' required.")
  }

  n_reps <- length(seeds)

  # Set up cluster
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  cat(sprintf("Running %d replications on %d cores...\n", n_reps, n_cores))

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))

  # Export functions and variables to workers
  parallel::clusterExport(cl, c(
    "run_one_rep",
    "ampute_data",
    "ampute_mcar_monotone",
    "ampute_mar_monotone",
    "ampute_mnar_monotone",
    "negloglik_complete",
    "negloglik_observed",
    "compute_Q_analytical",
    "compute_Q_analytical_proper",
    "estimate_imputation_params",
    "estimate_V_theta",
    "compute_hessian_correction_A",
    "compute_hessian_correction_B",
    "compute_hessian_correction_C",
    "seeds",
    "n",
    "mechanism",
    "gamma",
    "cut1",
    "cut2",
    "true_params",
    "imputation_model"
  ), envir = environment())

  # Run in parallel with progress bar
  results <- pbapply::pblapply(1:n_reps, function(rep) {
    run_one_rep(
      seed_data = seeds[[rep]]$data,
      seed_ampute = seeds[[rep]]$ampute,
      n = n,
      mechanism = mechanism,
      gamma = gamma,
      cut1 = cut1,
      cut2 = cut2,
      true_params = true_params,
      imputation_model = imputation_model
    )
  }, cl = cl)

  results <- do.call(rbind, Filter(Negate(is.null), results))
  as.data.frame(results)
}
