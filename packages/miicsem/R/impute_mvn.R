#' Conditional MVN Imputation Helper Given Implied (mu, Sigma)
#'
#' Workhorse for any model-based MVN imputation: given a target mean
#' vector and covariance matrix, draws Y_mis | Y_obs from the implied
#' multivariate normal conditional for every missing pattern in
#' \code{data_miss}, returning a list of M completed data frames.
#'
#' @param data_miss Amputed data.frame.
#' @param mu_hat Named numeric vector of length p.
#' @param Sigma_hat p x p covariance matrix (named rows/cols).
#' @param M Number of imputations.
#' @param seed Integer seed.
#' @param var_names Character vector ordering Y columns.
#' @return List of M completed data frames.
#' @keywords internal
impute_mvn_conditional <- function(data_miss, mu_hat, Sigma_hat,
                                   M, seed, var_names) {
  p <- length(var_names)

  # Reorder to var_names
  mu_hat    <- mu_hat[var_names]
  Sigma_hat <- Sigma_hat[var_names, var_names, drop = FALSE]

  # Identify unique missing patterns (rows).  Use a compact 0/1 string
  # so each position is one character, then parse back.
  miss_mat   <- is.na(as.matrix(data_miss[, var_names, drop = FALSE]))
  miss_01    <- matrix(as.integer(miss_mat), nrow = nrow(miss_mat))
  pattern_id <- apply(miss_01, 1, paste, collapse = "")
  unique_pat <- unique(pattern_id)

  # Precompute conditional distribution per unique missing pattern
  cond_spec <- lapply(unique_pat, function(pat) {
    pat_vec  <- as.integer(strsplit(pat, "")[[1]])
    miss_idx <- which(pat_vec == 1L)
    obs_idx  <- setdiff(seq_len(p), miss_idx)
    if (length(miss_idx) == 0) {
      return(list(miss_idx = integer(0), obs_idx = obs_idx,
                  mu_mis = NULL, b_mo = NULL, Sigma_cond = NULL,
                  Sigma_cond_chol = NULL, mu_obs = NULL))
    }
    if (length(obs_idx) == 0) {
      return(list(
        miss_idx = miss_idx, obs_idx = obs_idx,
        mu_mis   = mu_hat[miss_idx], b_mo = NULL,
        Sigma_cond = Sigma_hat[miss_idx, miss_idx, drop = FALSE],
        Sigma_cond_chol = tryCatch(
          chol(Sigma_hat[miss_idx, miss_idx, drop = FALSE]),
          error = function(e) NULL),
        mu_obs   = NULL
      ))
    }
    Sigma_oo <- Sigma_hat[obs_idx,  obs_idx,  drop = FALSE]
    Sigma_mo <- Sigma_hat[miss_idx, obs_idx,  drop = FALSE]
    Sigma_mm <- Sigma_hat[miss_idx, miss_idx, drop = FALSE]
    Sigma_oo_inv <- tryCatch(solve(Sigma_oo),
                             error = function(e) MASS::ginv(Sigma_oo))
    b_mo         <- Sigma_mo %*% Sigma_oo_inv
    Sigma_cond   <- Sigma_mm - b_mo %*% t(Sigma_mo)
    Sigma_cond <- (Sigma_cond + t(Sigma_cond)) / 2
    Sigma_cond_chol <- tryCatch(chol(Sigma_cond),
                                error = function(e) {
                                  d <- nrow(Sigma_cond)
                                  eps <- 1e-8 * mean(diag(Sigma_cond))
                                  chol(Sigma_cond + diag(eps, d))
                                })
    list(miss_idx = miss_idx, obs_idx = obs_idx,
         mu_mis = mu_hat[miss_idx], mu_obs = mu_hat[obs_idx],
         b_mo = b_mo, Sigma_cond = Sigma_cond,
         Sigma_cond_chol = Sigma_cond_chol)
  })
  names(cond_spec) <- unique_pat

  set.seed(seed)
  data_mat <- as.matrix(data_miss[, var_names, drop = FALSE])

  lapply(seq_len(M), function(m) {
    out <- data_mat
    for (pat in unique_pat) {
      spec <- cond_spec[[pat]]
      if (length(spec$miss_idx) == 0) next
      rows <- which(pattern_id == pat)
      if (is.null(spec$Sigma_cond_chol)) next
      if (length(spec$obs_idx) == 0) {
        Z <- matrix(stats::rnorm(length(rows) * length(spec$miss_idx)),
                    nrow = length(rows), ncol = length(spec$miss_idx))
        draws <- sweep(Z %*% spec$Sigma_cond_chol, 2, spec$mu_mis, "+")
      } else {
        y_obs   <- data_mat[rows, spec$obs_idx, drop = FALSE]
        centred <- sweep(y_obs, 2, spec$mu_obs, "-")
        cond_means <- sweep(centred %*% t(spec$b_mo), 2, spec$mu_mis, "+")
        Z <- matrix(stats::rnorm(length(rows) * length(spec$miss_idx)),
                    nrow = length(rows), ncol = length(spec$miss_idx))
        draws <- cond_means + Z %*% spec$Sigma_cond_chol
      }
      out[rows, spec$miss_idx] <- draws
    }
    df <- as.data.frame(out)
    colnames(df) <- var_names
    df
  })
}


#' Congenial Model-Based Imputation via Conditional Multivariate Normal
#'
#' Fits the saturated SEM model via FIML, then draws Y_mis | Y_obs from
#' the implied multivariate normal conditional via
#' \code{\link{impute_mvn_conditional}}.  Joint-coherent imputation,
#' MLE-plug-in (improper).
#'
#' @inheritParams impute_mvn_conditional
#' @return A list of \code{M} completed data.frames, or NULL if the
#'   FIML fit fails.
#' @export
impute_from_saturated_mvn <- function(data_miss, M, seed,
                                      var_names = NULL) {
  if (is.null(var_names)) var_names <- colnames(data_miss)
  sat_syntax <- get_saturated_model(var_names)
  fit_sat <- tryCatch({
    lavaan::sem(
      model         = sat_syntax,
      data          = data_miss,
      missing       = "fiml",
      meanstructure = TRUE,
      estimator     = "ML",
      se            = "none",
      test          = "standard",
      warn          = FALSE
    )
  }, error = function(e) NULL)
  if (is.null(fit_sat) || !lavaan::lavInspect(fit_sat, "converged")) {
    return(NULL)
  }

  mu_hat    <- lavaan::lavInspect(fit_sat, "mean.ov")
  Sigma_hat <- lavaan::lavInspect(fit_sat, "cov.ov")
  impute_mvn_conditional(data_miss, mu_hat, Sigma_hat,
                         M = M, seed = seed, var_names = var_names)
}


#' True-Model (M1) Imputation via Conditional Multivariate Normal
#'
#' Fits M1 (the true 3-factor CFA) via FIML on the amputed data, then
#' draws Y_mis | Y_obs from M1's *model-implied* (mu, Sigma).  This is
#' the strictest possible congenial imputation: imputation parameter
#' space = analysis parameter space (32 free parameters with the factor
#' structure built in), removing the 13 extra free directions in the
#' saturated MVN's Sigma.
#'
#' Useful as a theoretical-bound benchmark: under v4's congeniality
#' assumption with no extra slack, Term 1 should hit +tr(RIV) to leading
#' order.  Any residual gap from +1 reflects finite-N, finite-M, or
#' MLE-plug-in vs proper-MI effects, NOT imputation-model misspecification.
#'
#' @inheritParams impute_mvn_conditional
#' @return A list of \code{M} completed data.frames, or NULL if the
#'   M1 FIML fit fails.
#' @export
impute_from_M1_mvn <- function(data_miss, M, seed, var_names = NULL) {
  if (is.null(var_names)) var_names <- colnames(data_miss)
  M1_syntax <- get_sim1_models()$M1
  fit_M1 <- tryCatch({
    lavaan::sem(
      model         = M1_syntax,
      data          = data_miss,
      missing       = "fiml",
      meanstructure = TRUE,
      estimator     = "ML",
      se            = "none",
      test          = "standard",
      warn          = FALSE
    )
  }, error = function(e) NULL)
  if (is.null(fit_M1) || !lavaan::lavInspect(fit_M1, "converged")) {
    return(NULL)
  }

  # Model-implied moments under fitted M1 (NOT sample moments)
  implied   <- lavaan::fitted(fit_M1)
  mu_hat    <- implied$mean
  Sigma_hat <- implied$cov
  impute_mvn_conditional(data_miss, mu_hat, Sigma_hat,
                         M = M, seed = seed, var_names = var_names)
}
