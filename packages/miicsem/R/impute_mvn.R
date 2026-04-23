#' Congenial Model-Based Imputation via Conditional Multivariate Normal
#'
#' For each missing pattern, draws Y_mis | Y_obs from the multivariate
#' normal conditional implied by a fit of the saturated model \code{Msat}
#' to the amputed data via FIML.  Because the saturated model is
#' congenial with any CFA candidate nested within it (all share the same
#' multivariate normal family), this is a fully congenial imputation
#' strategy under the null that the data are multivariate normal.
#'
#' Mechanics: fit Msat via FIML -> obtain (mu_hat, Sigma_hat).  For each
#' unique missing pattern, partition into observed/missing blocks and
#' apply standard Gaussian conditional formulas:
#'
#'   mu_{mis|obs}(y_obs) = mu_mis + Sigma_{mo} Sigma_{oo}^{-1} (y_obs - mu_obs)
#'   Sigma_{mis|obs}     = Sigma_mm - Sigma_{mo} Sigma_{oo}^{-1} Sigma_{om}
#'
#' then draw Y_mis ~ N(mu_{mis|obs}(y_obs), Sigma_{mis|obs}) for each case
#' with that pattern, independently across M imputations.
#'
#' @param data_miss Amputed data.frame (may contain NAs).
#' @param M Number of imputations.
#' @param seed Integer seed for reproducibility.
#' @param var_names Character vector of variable names (for sanity
#'   checks; inferred from \code{data_miss} if NULL).
#' @return A list of \code{M} completed data.frames with the same column
#'   order and row order as \code{data_miss}.  Returns \code{NULL} on
#'   fatal failure (e.g. Msat FIML fit fails).
#' @export
impute_from_saturated_mvn <- function(data_miss, M, seed,
                                      var_names = NULL) {
  if (is.null(var_names)) var_names <- colnames(data_miss)
  p <- length(var_names)

  # Fit saturated model via FIML to get mu_hat, Sigma_hat
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
  # Ensure the orderings match var_names
  mu_hat    <- mu_hat[var_names]
  Sigma_hat <- Sigma_hat[var_names, var_names, drop = FALSE]

  # Identify unique missing patterns (rows).  Use a compact 0/1 string so
  # each position is one character, then parse back to a logical vector.
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
      # All missing: use marginal
      return(list(
        miss_idx = miss_idx,
        obs_idx  = obs_idx,
        mu_mis   = mu_hat[miss_idx],
        b_mo     = NULL,
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
    # Ensure symmetry (numerical)
    Sigma_cond <- (Sigma_cond + t(Sigma_cond)) / 2
    Sigma_cond_chol <- tryCatch(chol(Sigma_cond),
                                error = function(e) {
                                  # Add tiny ridge, retry
                                  d <- nrow(Sigma_cond)
                                  eps <- 1e-8 * mean(diag(Sigma_cond))
                                  chol(Sigma_cond + diag(eps, d))
                                })
    list(miss_idx   = miss_idx,
         obs_idx    = obs_idx,
         mu_mis     = mu_hat[miss_idx],
         mu_obs     = mu_hat[obs_idx],
         b_mo       = b_mo,
         Sigma_cond = Sigma_cond,
         Sigma_cond_chol = Sigma_cond_chol)
  })
  names(cond_spec) <- unique_pat

  # Generate M imputations
  set.seed(seed)
  data_mat <- as.matrix(data_miss[, var_names, drop = FALSE])

  imputed_list <- lapply(seq_len(M), function(m) {
    out <- data_mat
    for (pat in unique_pat) {
      spec <- cond_spec[[pat]]
      if (length(spec$miss_idx) == 0) next
      rows <- which(pattern_id == pat)
      if (is.null(spec$Sigma_cond_chol)) next

      if (length(spec$obs_idx) == 0) {
        # Marginal draws
        Z <- matrix(stats::rnorm(length(rows) * length(spec$miss_idx)),
                    nrow = length(rows),
                    ncol = length(spec$miss_idx))
        draws <- sweep(Z %*% spec$Sigma_cond_chol, 2,
                       spec$mu_mis, "+")
      } else {
        # Conditional draws per row
        y_obs  <- data_mat[rows, spec$obs_idx, drop = FALSE]
        centred <- sweep(y_obs, 2, spec$mu_obs, "-")
        cond_means <- sweep(centred %*% t(spec$b_mo), 2,
                            spec$mu_mis, "+")
        Z <- matrix(stats::rnorm(length(rows) * length(spec$miss_idx)),
                    nrow = length(rows),
                    ncol = length(spec$miss_idx))
        draws <- cond_means + Z %*% spec$Sigma_cond_chol
      }
      out[rows, spec$miss_idx] <- draws
    }
    df <- as.data.frame(out)
    colnames(df) <- var_names
    df
  })

  imputed_list
}
