# ============================================================================
# Utility Functions for 00-premise-basic-Y.R
# ============================================================================
# Functions for:
#   - Grid search for selection entropy calibration
#   - Imputation bias validation
#   - Visualization of results
# ============================================================================
# NOTE: Parallel processing should be configured in the calling script
#       before calling these functions (using setup_parallel())
# ============================================================================

#' Calibrate Selection Entropy
#'
#' Performs grid search over (n, rho) to find settings with good model
#' selection entropy (P(correct) ≈ 0.5-0.7, entropy ≈ 0.9-1.1)
#'
#' @param p Integer. Dimension
#' @param true_structure Character. True covariance structure
#' @param n_grid Integer vector. Sample sizes to test
#' @param rho_grid Numeric vector. Correlation parameters to test
#' @param candidate_models Character vector. Models to compare
#' @param n_reps Integer. Number of Monte Carlo replications per grid point
#' @param base_seed Integer. Base random seed
#' @param use_parallel Logical. Use parallel processing (assumes already configured)
#' @return data.frame with grid results
calibrate_selection_entropy <- function(
  p = 10,
  true_structure = "Toeplitz",
  n_grid = c(50, 100, 200, 500),
  rho_grid = c(0.3, 0.5, 0.7),
  candidate_models = c("CS", "Toeplitz", "Unstructured"),
  n_reps = 1000,
  base_seed = 12345,
  use_parallel = TRUE
) {

  cat("=== Grid Search for Selection Entropy Calibration ===\n")
  cat(sprintf("True structure: %s\n", true_structure))
  cat(sprintf("Dimension: p = %d\n", p))
  cat(sprintf("Grid: n ∈ {%s}, rho ∈ {%s}\n",
              paste(n_grid, collapse = ", "),
              paste(rho_grid, collapse = ", ")))
  cat(sprintf("Replications: %d\n", n_reps))
  cat(sprintf("Parallel: %s\n\n", ifelse(use_parallel, "Yes", "No")))

  # Create grid of (n, rho) combinations
  grid_params <- expand.grid(
    n = n_grid,
    rho = rho_grid,
    stringsAsFactors = FALSE
  )
  n_grid_points <- nrow(grid_params)

  cat(sprintf("Total grid points: %d\n", n_grid_points))
  cat(sprintf("Total datasets to generate: %d\n\n", n_grid_points * n_reps))

  # Function to run one grid point
  run_grid_point <- function(grid_idx, seed_offset) {

    n_val <- grid_params$n[grid_idx]
    rho_val <- grid_params$rho[grid_idx]

    # Generate true covariance
    Sigma_true <- generate_covariance(
      p = p,
      structure = true_structure,
      sigma2 = 1,
      rho = rho_val,
      seed = base_seed + grid_idx * 1000
    )
    mu_true <- rep(0, p)

    # Generate seeds for replications
    rep_seeds <- generate_seeds(n_reps, base_seed = seed_offset + grid_idx * 10000)

    # Run replications for this grid point
    selections <- sapply(1:n_reps, function(r) {

      # Generate complete data
      data_complete <- generate_mvn_data(
        n = n_val,
        mu = mu_true,
        Sigma = Sigma_true,
        seed = rep_seeds[r]
      )

      # Fit all candidate models and compute BIC
      bic_vals <- sapply(candidate_models, function(model) {
        fit_result <- fit_lavaan_model(data_complete, structure = model, return_vcov = FALSE)
        if (fit_result$converged) {
          return(fit_result$BIC)
        } else {
          return(Inf)  # Penalize non-convergence
        }
      })

      # Select model with minimum BIC
      selected_model <- candidate_models[which.min(bic_vals)]
      return(selected_model)
    })

    # Compute selection frequencies and entropy
    freq_summary <- compute_selection_frequencies(selections, true_model = true_structure)

    result <- data.frame(
      n = n_val,
      rho = rho_val,
      true_structure = true_structure,
      p_correct = attr(freq_summary, "p_correct"),
      entropy = attr(freq_summary, "entropy"),
      n_reps = n_reps
    )

    # Add individual model selection proportions
    for (model in candidate_models) {
      prop_col <- paste0("prop_", model)
      if (model %in% freq_summary$model) {
        result[[prop_col]] <- freq_summary$proportion[freq_summary$model == model]
      } else {
        result[[prop_col]] <- 0
      }
    }

    return(result)
  }

  # Run grid search (parallel or sequential based on use_parallel flag)
  if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    cat("Running grid search in parallel...\n")
    grid_results_list <- future.apply::future_lapply(1:n_grid_points, function(i) {
      run_grid_point(i, seed_offset = base_seed)
    }, future.seed = TRUE)
  } else {
    if (use_parallel) {
      warning("Parallel processing requested but not available; running sequentially")
    }
    cat("Running grid search sequentially...\n")
    grid_results_list <- lapply(1:n_grid_points, function(i) {
      run_grid_point(i, seed_offset = base_seed)
    })
  }

  # Combine results
  grid_results <- do.call(rbind, grid_results_list)

  cat("\n=== Grid Search Complete ===\n")
  cat("Summary of results:\n")
  print(grid_results[, c("n", "rho", "p_correct", "entropy")])

  # Identify good calibration points
  good_points <- grid_results %>%
    dplyr::filter(p_correct >= 0.5 & p_correct <= 0.7,
           entropy >= 0.9 & entropy <= 1.1)

  if (nrow(good_points) > 0) {
    cat("\n\nRecommended calibration points (P(correct) ≈ 0.6, entropy ≈ 1.0):\n")
    print(good_points[, c("n", "rho", "p_correct", "entropy")])
  } else {
    cat("\n\nNo grid points met target criteria. Closest:\n")
    grid_results$dist_target <- abs(grid_results$p_correct - 0.6) + abs(grid_results$entropy - 1.0)
    closest <- grid_results[order(grid_results$dist_target)[1:min(3, nrow(grid_results))], ]
    print(closest[, c("n", "rho", "p_correct", "entropy")])
  }

  return(grid_results)
}


#' Validate Imputation Bias
#'
#' Empirically validates that imputation bias ≈ tr(RIV)
#'
#' @param n Integer. Sample size
#' @param p Integer. Dimension
#' @param true_structure Character. True covariance structure
#' @param rho Numeric. Correlation parameter
#' @param missing_rate Numeric. Proportion of values missing
#' @param M Integer. Number of imputations
#' @param scenarios Character vector. Imputation scenarios ("True", "MLE")
#' @param n_reps Integer. Number of Monte Carlo replications
#' @param base_seed Integer. Base random seed
#' @param use_parallel Logical. Use parallel processing (assumes already configured)
#' @return data.frame with bias validation results
validate_imputation_bias <- function(
  n = 200,
  p = 10,
  true_structure = "Toeplitz",
  rho = 0.5,
  missing_rate = 0.6,
  M = 100,
  scenarios = c("True", "MLE"),
  n_reps = 1000,
  base_seed = 54321,
  use_parallel = TRUE
) {

  cat("\n=== Imputation Bias Validation ===\n")
  cat(sprintf("True structure: %s (p = %d)\n", true_structure, p))
  cat(sprintf("Sample size: n = %d\n", n))
  cat(sprintf("Missing rate: %.1f%%\n", missing_rate * 100))
  cat(sprintf("Imputations: M = %d\n", M))
  cat(sprintf("Scenarios: %s\n", paste(scenarios, collapse = ", ")))
  cat(sprintf("Replications: %d\n\n", n_reps))

  # Generate true parameters
  Sigma_true <- generate_covariance(
    p = p,
    structure = true_structure,
    sigma2 = 1,
    rho = rho,
    seed = base_seed
  )
  mu_true <- rep(0, p)

  # Get true Q
  Q_true <- get_param_count(p = p, structure = true_structure)
  cat(sprintf("True model complexity: Q = %d\n\n", Q_true))

  # Function to run one replication
  run_replication <- function(rep_idx, seed_val) {

    # Generate complete data
    data_complete <- generate_mvn_data(
      n = n,
      mu = mu_true,
      Sigma = Sigma_true,
      seed = seed_val
    )

    # Impose missingness
    miss_result <- impose_missingness(
      data = data_complete,
      missing_rate = missing_rate,
      pattern = "monotone",
      prop_complete = 0.4,
      seed = seed_val + 1
    )
    data_miss <- miss_result$data_miss

    # Initialize result storage for this replication
    rep_results <- list()

    # Loop over scenarios
    for (scenario in scenarios) {

      # Generate M imputations
      if (scenario == "True") {
        impute_result <- impute_wrapper(
          data_miss = data_miss,
          scenario = "True",
          mu_true = mu_true,
          Sigma_true = Sigma_true,
          M = M,
          seed_base = seed_val + 100
        )
      } else if (scenario == "MLE") {
        impute_result <- impute_wrapper(
          data_miss = data_miss,
          scenario = "MLE",
          structure = true_structure,
          M = M,
          seed_base = seed_val + 200
        )
      }

      completed_datasets <- impute_result$completed_datasets

      # Fit model to observed data to get θ̂_obs
      fit_obs <- fit_lavaan_model(data_miss, structure = true_structure, return_vcov = TRUE)
      if (!fit_obs$converged) {
        warning(paste("Replication", rep_idx, "scenario", scenario, ": fit did not converge"))
      }

      # Compute Rubin's variance components and RIV
      riv_result <- tryCatch(
        {
          compute_rubin_variance(completed_datasets, structure = true_structure)
        },
        error = function(e) {
          warning(paste("Replication", rep_idx, "scenario", scenario,
                       ": compute_rubin_variance failed -", e$message))
          return(NULL)
        }
      )

      # Skip this scenario if RIV computation failed
      if (is.null(riv_result)) {
        next
      }

      tr_RIV <- riv_result$tr_RIV
      M_successful <- riv_result$M_successful

      # Compute Q̄_MI(θ̂_obs)
      theta_hat_list <- list(mu_hat = fit_obs$mu_hat, Sigma_hat = fit_obs$Sigma_hat)
      Q_MI <- compute_Q_function(theta_hat_list, completed_datasets)

      # Compute ℓ_com(θ̂_obs) on true complete data
      ell_com <- compute_complete_loglik(data_complete, fit_obs$mu_hat, fit_obs$Sigma_hat)

      # Compute imputation bias
      bias_empirical <- compute_imputation_bias(Q_MI, ell_com)

      # Store results
      rep_results[[scenario]] <- data.frame(
        replication = rep_idx,
        scenario = scenario,
        n = n,
        p = p,
        Q = Q_true,
        M = M,
        M_successful = M_successful,
        missing_rate = missing_rate,
        true_structure = true_structure,
        bias_empirical = bias_empirical,
        tr_RIV_theoretical = tr_RIV,
        Q_MI = Q_MI,
        ell_com = ell_com,
        converged = fit_obs$converged,
        stringsAsFactors = FALSE
      )
    }

    # Combine scenarios for this replication
    rep_df <- do.call(rbind, rep_results)
    return(rep_df)
  }

  # Generate seeds for replications
  rep_seeds <- generate_seeds(n_reps, base_seed = base_seed)

  # Run replications (parallel or sequential based on use_parallel flag)
  if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    cat("Running validation replications in parallel...\n")
    results_list <- future.apply::future_lapply(1:n_reps, function(i) {
      run_replication(i, rep_seeds[i])
    }, future.seed = TRUE)
  } else {
    if (use_parallel) {
      warning("Parallel processing requested but not available; running sequentially")
    }
    cat("Running validation replications sequentially...\n")
    results_list <- lapply(1:n_reps, function(i) {
      run_replication(i, rep_seeds[i])
    })
  }

  # Combine results
  validation_df <- do.call(rbind, results_list)

  # Print summary
  cat("\n=== Validation Complete ===\n")
  summary_by_scenario <- validation_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      mean_M_successful = mean(M_successful, na.rm = TRUE),
      min_M_successful = min(M_successful, na.rm = TRUE),
      mean_bias = mean(bias_empirical, na.rm = TRUE),
      sd_bias = sd(bias_empirical, na.rm = TRUE),
      mean_tr_RIV = mean(tr_RIV_theoretical, na.rm = TRUE),
      correlation = cor(bias_empirical, tr_RIV_theoretical, use = "complete.obs"),
      .groups = "drop"
    )

  print(summary_by_scenario)

  # Report imputation success rate
  overall_success_rate <- mean(validation_df$M_successful / validation_df$M)
  cat(sprintf("\nOverall imputation success rate: %.1f%%\n",
              overall_success_rate * 100))

  return(validation_df)
}


#' Visualize Grid Search Results
#'
#' Creates heatmaps for P(correct) and entropy by (n, rho)
#'
#' @param grid_results data.frame. Output from calibrate_selection_entropy()
#' @param output_dir Character. Directory for saving figures
visualize_grid_search <- function(grid_results, output_dir = "simulations/00-premise-basic-Y/figures") {

  cat("\n=== Generating Grid Search Visualizations ===\n")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Heatmap: P(correct model selected)
  p1 <- ggplot2::ggplot(grid_results, ggplot2::aes(x = factor(n), y = factor(rho), fill = p_correct)) +
    ggplot2::geom_tile(color = "white", size = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", p_correct)), color = "black", size = 4) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                         midpoint = 0.6, limits = c(0, 1)) +
    ggplot2::labs(title = "Model Selection Accuracy by (n, ρ)",
         subtitle = sprintf("True model: %s", unique(grid_results$true_structure)[1]),
         x = "Sample Size (n)", y = "Correlation (ρ)",
         fill = "P(Correct)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(file.path(output_dir, "grid_heatmap_accuracy.pdf"), p1, width = 8, height = 6)
  ggplot2::ggsave(file.path(output_dir, "grid_heatmap_accuracy.png"), p1, width = 8, height = 6)

  # Heatmap: Entropy
  p2 <- ggplot2::ggplot(grid_results, ggplot2::aes(x = factor(n), y = factor(rho), fill = entropy)) +
    ggplot2::geom_tile(color = "white", size = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", entropy)), color = "black", size = 4) +
    ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "blue",
                         midpoint = 1.0, limits = c(0, max(grid_results$entropy))) +
    ggplot2::labs(title = "Selection Entropy by (n, ρ)",
         subtitle = "Target: entropy ≈ 0.9-1.1 for interesting selection problem",
         x = "Sample Size (n)", y = "Correlation (ρ)",
         fill = "Entropy") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(file.path(output_dir, "grid_heatmap_entropy.pdf"), p2, width = 8, height = 6)
  ggplot2::ggsave(file.path(output_dir, "grid_heatmap_entropy.png"), p2, width = 8, height = 6)

  cat(sprintf("Saved grid search visualizations to %s\n", output_dir))
}


#' Visualize Bias Validation Results
#'
#' Creates diagnostic plots for bias validation
#'
#' @param validation_df data.frame. Output from validate_imputation_bias()
#' @param output_dir Character. Directory for saving figures
visualize_bias_validation <- function(validation_df, output_dir = "simulations/00-premise-basic-Y/figures") {

  cat("\n=== Generating Bias Validation Visualizations ===\n")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 1. Scatterplot: Bias_empirical vs tr(RIV)_theoretical
  p1 <- ggplot2::ggplot(validation_df, ggplot2::aes(x = tr_RIV_theoretical, y = bias_empirical, color = scenario)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
    ggplot2::facet_wrap(~ true_structure + scenario, scales = "free") +
    ggplot2::labs(title = "Empirical Bias vs Theoretical tr(RIV)",
         subtitle = "45° line: perfect agreement",
         x = "tr(RIV) [Theoretical]",
         y = "Bias [Empirical]",
         color = "Scenario") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "validation_bias_vs_RIV.pdf"), p1, width = 10, height = 6)
  ggplot2::ggsave(file.path(output_dir, "validation_bias_vs_RIV.png"), p1, width = 10, height = 6)

  # 2. Boxplot: Bias by scenario
  p2 <- ggplot2::ggplot(validation_df, ggplot2::aes(x = scenario, y = bias_empirical, fill = scenario)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~ true_structure) +
    ggplot2::labs(title = "Imputation Bias by Scenario",
         subtitle = "True scenario should have bias ≈ 0",
         x = "Scenario", y = "Imputation Bias",
         fill = "Scenario") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "validation_bias_by_scenario.pdf"), p2, width = 10, height = 6)
  ggplot2::ggsave(file.path(output_dir, "validation_bias_by_scenario.png"), p2, width = 10, height = 6)

  # 3. Scatterplot: Bias vs Q
  p3 <- ggplot2::ggplot(validation_df, ggplot2::aes(x = Q, y = bias_empirical, color = true_structure)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_smooth(method = "lm", se = TRUE) +
    ggplot2::facet_wrap(~ scenario) +
    ggplot2::labs(title = "Bias Scaling with Model Complexity",
         subtitle = "Bias should increase linearly with Q",
         x = "Model Complexity (Q)",
         y = "Imputation Bias",
         color = "Structure") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "validation_bias_vs_Q.pdf"), p3, width = 10, height = 6)
  ggplot2::ggsave(file.path(output_dir, "validation_bias_vs_Q.png"), p3, width = 10, height = 6)

  # 4. Export summary table
  summary_table <- validation_df %>%
    dplyr::group_by(true_structure, p, Q, M, scenario) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      mean_bias = mean(bias_empirical, na.rm = TRUE),
      sd_bias = sd(bias_empirical, na.rm = TRUE),
      mean_tr_RIV = mean(tr_RIV_theoretical, na.rm = TRUE),
      correlation = cor(bias_empirical, tr_RIV_theoretical, use = "complete.obs"),
      .groups = "drop"
    )

  utils::write.csv(summary_table,
            file.path(output_dir, "validation_summary_table.csv"),
            row.names = FALSE)

  cat(sprintf("Saved bias validation visualizations to %s\n", output_dir))
}


# ============================================================================
# Term Decomposition Validation
# ============================================================================
# Empirically verify Terms 1-3 from the Taylor expansion in the derivation
# to identify the source of the 2×tr(RIV) discrepancy.
# ============================================================================

#' Validate Term Decomposition of Imputation Bias
#'
#' Empirically computes Terms 1-3 from the Taylor expansion to verify which
#' term(s) contribute to the observed 2×tr(RIV) bias.
#'
#' @param n Integer. Sample size
#' @param p Integer. Number of variables
#' @param true_structure Character. "CS", "Toeplitz", or "Unstructured"
#' @param rho Numeric. Correlation parameter
#' @param missing_rate Numeric. Proportion of missing data
#' @param M Integer. Number of imputations
#' @param n_reps Integer. Number of replications
#' @param base_seed Integer. Base seed for reproducibility
#' @param use_parallel Logical. Whether to use parallel processing
#' @return data.frame with columns for each term and diagnostic values
#' @details
#'   For each replication, computes:
#'   - Total = Q̄_MI(θ̂) - ℓ_com(θ̂) [the bias we observe]
#'   - Term1 = Q̄_MI(θ*) - ℓ_com(θ*) [should be ≈ 0]
#'   - Term2 = (θ̂-θ*)ᵀ[score_diff] [should be ≈ tr(RIV)]
#'   - Term3 = -½(θ̂-θ*)ᵀ[hessian_diff](θ̂-θ*) [should be negligible]
#'   - Taylor Sum = Term1 + Term2 + Term3
#'   - Residual = Total - Taylor Sum
#'
#'   If Taylor Sum ≈ Total, the expansion is valid. Then examine which term(s)
#'   account for the extra tr(RIV) observed in the bias.
validate_term_decomposition <- function(
  n = 200,
  p = 10,
  true_structure = "Toeplitz",
  rho = 0.5,
  missing_rate = 0.6,
  M = 100,
  n_reps = 100,
  base_seed = 54321,
  use_parallel = FALSE
) {

  cat("=== Term Decomposition Validation ===\n")
  cat(sprintf("n=%d, p=%d, structure=%s, rho=%.2f, missing=%.0f%%, M=%d, reps=%d\n",
              n, p, true_structure, rho, missing_rate * 100, M, n_reps))

  # Generate true parameters
  Sigma_true <- generate_covariance(p = p, structure = true_structure, sigma2 = 1, rho = rho)
  mu_true <- rep(0, p)

  # Compute Q (number of parameters) for the unstructured natural parameterization
  Q_natural <- p + p * (p + 1) / 2  # μ + vech(Σ)

  cat(sprintf("True parameters: Q_natural=%d (μ=%d, Σ=%d)\n",
              as.integer(Q_natural), p, as.integer(p * (p + 1) / 2)))

  # Inner function for one replication
  run_replication <- function(rep_idx, seed_val) {
    # Set seed for this replication
    set.seed(seed_val)

    # Generate complete data
    data_complete <- generate_mvn_data(n, mu_true, Sigma_true, seed = seed_val)

    # Impose missingness
    miss_result <- impose_missingness(
      data = data_complete,
      missing_rate = missing_rate,
      pattern = "monotone",
      prop_complete = 0.4,
      seed = seed_val + 1000
    )
    data_miss <- miss_result$data_miss

    # Generate imputations using MLE parameters
    impute_result <- impute_wrapper(
      data_miss = data_miss,
      scenario = "MLE",
      structure = true_structure,
      M = M,
      seed_base = seed_val + 2000
    )
    completed_datasets <- impute_result$completed_datasets

    # Fit model to observed data to get θ̂_obs
    fit_obs <- fit_lavaan_model(data_miss, structure = true_structure, return_vcov = TRUE)
    if (!fit_obs$converged) {
      warning(paste("Replication", rep_idx, ": observed-data fit did not converge"))
      return(NULL)
    }

    mu_hat <- fit_obs$mu_hat
    Sigma_hat <- fit_obs$Sigma_hat

    # Compute tr(RIV) via Rubin's variance components
    riv_result <- tryCatch(
      {
        compute_rubin_variance(completed_datasets, structure = true_structure)
      },
      error = function(e) {
        warning(paste("Replication", rep_idx, ": compute_rubin_variance failed -", e$message))
        return(NULL)
      }
    )

    if (is.null(riv_result)) {
      return(NULL)
    }

    tr_RIV <- riv_result$tr_RIV

    # Compute TOTAL bias: Q̄_MI(θ̂) - ℓ_com(θ̂)
    theta_hat_list <- list(mu_hat = mu_hat, Sigma_hat = Sigma_hat)
    Q_MI_at_hat <- compute_Q_function(theta_hat_list, completed_datasets)
    ell_com_at_hat <- compute_complete_loglik(data_complete, mu_hat, Sigma_hat)
    total_bias <- Q_MI_at_hat - ell_com_at_hat

    # Compute Term 1: Q̄_MI(θ*) - ℓ_com(θ*)
    term1 <- tryCatch(
      {
        compute_term1(completed_datasets, data_complete, mu_true, Sigma_true)
      },
      error = function(e) {
        warning(paste("Replication", rep_idx, ": compute_term1 failed -", e$message))
        return(NA)
      }
    )

    # Compute Term 2: (θ̂ - θ*)ᵀ [E[S_com(θ*)|θ̃] - S_com(θ*)]
    term2 <- tryCatch(
      {
        compute_term2(completed_datasets, data_complete, mu_true, Sigma_true, mu_hat, Sigma_hat)
      },
      error = function(e) {
        warning(paste("Replication", rep_idx, ": compute_term2 failed -", e$message))
        return(NA)
      }
    )

    # Compute Term 3: -½(θ̂ - θ*)ᵀ [E[J_com(θ*)|θ̃] - J_com(θ*)] (θ̂ - θ*)
    term3 <- tryCatch(
      {
        compute_term3(completed_datasets, data_complete, mu_true, Sigma_true, mu_hat, Sigma_hat)
      },
      error = function(e) {
        warning(paste("Replication", rep_idx, ": compute_term3 failed -", e$message))
        return(NA)
      }
    )

    # Taylor sum and residual
    taylor_sum <- term1 + term2 + term3
    residual <- total_bias - taylor_sum

    return(data.frame(
      replication = rep_idx,
      n = n,
      p = p,
      Q_natural = Q_natural,
      M = M,
      missing_rate = missing_rate,
      true_structure = true_structure,
      total_bias = total_bias,
      term1 = term1,
      term2 = term2,
      term3 = term3,
      taylor_sum = taylor_sum,
      residual = residual,
      tr_RIV = tr_RIV,
      converged = fit_obs$converged,
      stringsAsFactors = FALSE
    ))
  }

  # Generate seeds for replications
  rep_seeds <- generate_seeds(n_reps, base_seed = base_seed)

  # Run replications with progress bar
  if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    cat("Running term decomposition in parallel...\n")
    cat("(Progress bar not available in parallel mode)\n")

    results_list <- future.apply::future_lapply(1:n_reps, function(i) {
      run_replication(i, rep_seeds[i])
    }, future.seed = TRUE)
  } else {
    if (use_parallel) {
      warning("Parallel processing requested but not available; running sequentially")
    }
    cat("Running term decomposition sequentially...\n")

    # Use progressr for progress bar if available
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::handlers("txtprogressbar")

      results_list <- progressr::with_progress({
        p <- progressr::progressor(steps = n_reps)
        lapply(1:n_reps, function(i) {
          result <- run_replication(i, rep_seeds[i])
          p()  # Signal progress after completion
          result
        })
      })
    } else {
      results_list <- lapply(1:n_reps, function(i) {
        run_replication(i, rep_seeds[i])
      })
    }
  }

  # Filter out failed replications
  results_list <- results_list[!sapply(results_list, is.null)]
  n_success <- length(results_list)

  if (n_success == 0) {
    stop("All replications failed")
  }

  if (n_success < n_reps) {
    warning(sprintf("%d out of %d replications failed", n_reps - n_success, n_reps))
  }

  # Combine results
  validation_df <- do.call(rbind, results_list)

  # Print summary
  cat("\n=== Term Decomposition Summary ===\n")
  cat(sprintf("Successful replications: %d / %d\n\n", n_success, n_reps))

  # Compute summary statistics
  mean_tr_RIV <- mean(validation_df$tr_RIV, na.rm = TRUE)

  summary_stats <- data.frame(
    Metric = c("Total Bias", "Term 1", "Term 2", "Term 3", "Taylor Sum", "Residual", "tr(RIV)"),
    Mean = c(
      mean(validation_df$total_bias, na.rm = TRUE),
      mean(validation_df$term1, na.rm = TRUE),
      mean(validation_df$term2, na.rm = TRUE),
      mean(validation_df$term3, na.rm = TRUE),
      mean(validation_df$taylor_sum, na.rm = TRUE),
      mean(validation_df$residual, na.rm = TRUE),
      mean_tr_RIV
    ),
    SD = c(
      sd(validation_df$total_bias, na.rm = TRUE),
      sd(validation_df$term1, na.rm = TRUE),
      sd(validation_df$term2, na.rm = TRUE),
      sd(validation_df$term3, na.rm = TRUE),
      sd(validation_df$taylor_sum, na.rm = TRUE),
      sd(validation_df$residual, na.rm = TRUE),
      sd(validation_df$tr_RIV, na.rm = TRUE)
    )
  )

  # Add ratio to tr(RIV)
  summary_stats$Ratio_to_trRIV <- summary_stats$Mean / mean_tr_RIV

  print(summary_stats)

  cat("\n=== Interpretation ===\n")
  term2_ratio <- mean(validation_df$term2, na.rm = TRUE) / mean_tr_RIV
  term3_ratio <- mean(validation_df$term3, na.rm = TRUE) / mean_tr_RIV
  residual_ratio <- abs(mean(validation_df$residual, na.rm = TRUE)) / abs(mean(validation_df$total_bias, na.rm = TRUE))

  cat(sprintf("Term 2 / tr(RIV) = %.2f (theory predicts 1.0)\n", term2_ratio))
  cat(sprintf("Term 3 / tr(RIV) = %.2f (theory predicts ~0)\n", term3_ratio))
  cat(sprintf("Residual / Total = %.2f%% (Taylor approximation error)\n", residual_ratio * 100))

  if (abs(term3_ratio) > 0.5) {
    cat("\n** Term 3 is NOT negligible - Step 15 may be incorrect **\n")
  }

  if (term2_ratio > 1.5) {
    cat("\n** Term 2 > 1.5×tr(RIV) - Steps 14-16 may have an error **\n")
  }

  return(validation_df)
}
