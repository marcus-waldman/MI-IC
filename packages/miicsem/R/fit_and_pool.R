#' Fit All Candidate + Saturated Models to Complete Data (Oracle)
#'
#' @param data_complete Complete data.frame (no missingness).
#' @param models Named list of lavaan model syntaxes.
#' @param pop_starts Optional named list of warm-start coefficient
#'   vectors (as produced by \code{\link{compute_pop_starts}}). Passed
#'   via lavaan's \code{start} argument.
#' @return Named list; each element has \code{loglik}, \code{npar},
#'   \code{converged}.
fit_complete <- function(data_complete, models, pop_starts = NULL) {
  results <- lapply(names(models), function(mname) {
    tryCatch({
      start_val <- if (!is.null(pop_starts) && !is.null(pop_starts[[mname]])) {
        pop_starts[[mname]]
      } else {
        "default"
      }
      fit <- lavaan::sem(
        model         = models[[mname]],
        data          = data_complete,
        meanstructure = TRUE,
        estimator     = "ML",
        se            = "standard",
        test          = "standard",
        start         = start_val,
        warn          = FALSE
      )
      list(
        loglik            = as.numeric(lavaan::logLik(fit)),
        npar              = lavaan::fitmeasures(fit, "npar"),
        converged         = lavaan::lavInspect(fit, "converged"),
        vcov              = lavaan::vcov(fit),
        coef              = lavaan::coef(fit),
        partable_template = lavaan::parameterTable(fit)
      )
    }, error = function(e) {
      list(loglik = NA, npar = NA, converged = FALSE,
           vcov = NULL, coef = NULL, partable_template = NULL)
    })
  })
  names(results) <- names(models)
  return(results)
}


#' Fit All Candidate + Saturated Models on Amputed Data via FIML
#'
#' Provides the observed-data vcov V_obs at the FIML MLE, which is the
#' ingredient for the theoretical tr(RIV) reference:
#' \code{tr(V_obs \%*\% solve(V_com)) - p}.
#'
#' @param data_miss Amputed data.frame (pre-imputation).
#' @param models Named list of lavaan model syntaxes.
#' @param pop_starts Optional named list of warm-start coefficient
#'   vectors from \code{\link{compute_pop_starts}}.
#' @return Named list; each element has \code{loglik}, \code{npar},
#'   \code{converged}, \code{vcov}, \code{coef}.
fit_observed <- function(data_miss, models, pop_starts = NULL) {
  results <- lapply(names(models), function(mname) {
    tryCatch({
      start_val <- if (!is.null(pop_starts) && !is.null(pop_starts[[mname]])) {
        pop_starts[[mname]]
      } else {
        "default"
      }
      fit <- lavaan::sem(
        model         = models[[mname]],
        data          = data_miss,
        missing       = "fiml",
        meanstructure = TRUE,
        estimator     = "ML",
        se            = "standard",
        test          = "standard",
        start         = start_val,
        warn          = FALSE
      )
      list(
        loglik    = as.numeric(lavaan::logLik(fit)),
        npar      = lavaan::fitmeasures(fit, "npar"),
        converged = lavaan::lavInspect(fit, "converged"),
        vcov      = lavaan::vcov(fit),
        coef      = lavaan::coef(fit)
      )
    }, error = function(e) {
      list(loglik = NA, npar = NA, converged = FALSE,
           vcov = NULL, coef = NULL)
    })
  })
  names(results) <- names(models)
  return(results)
}


#' Fit a Single Model to All M Imputed Datasets
#'
#' @param imputed_list List of M data.frames.
#' @param model_syntax lavaan model syntax (character).
#' @param start_coefs Optional named numeric vector of starting
#'   coefficients (e.g., from \code{\link{compute_pop_starts}}). Passed
#'   via lavaan's \code{start} argument to every per-imputation fit.
#' @return List with \code{coefs}, \code{vcovs}, \code{logliks},
#'   \code{npar}, \code{param_names}, \code{converged_count},
#'   \code{converged_flags}, \code{partable_template}, \code{success}.
fit_single_model_mi <- function(imputed_list, model_syntax, start_coefs = NULL) {
  M <- length(imputed_list)
  start_val <- if (!is.null(start_coefs)) start_coefs else "default"

  fits <- lapply(seq_len(M), function(m) {
    tryCatch({
      fit <- lavaan::sem(
        model         = model_syntax,
        data          = imputed_list[[m]],
        meanstructure = TRUE,
        estimator     = "ML",
        se            = "standard",
        test          = "standard",
        start         = start_val,
        warn          = FALSE
      )
      list(
        coef      = lavaan::coef(fit),
        vcov      = lavaan::vcov(fit),
        loglik    = as.numeric(lavaan::logLik(fit)),
        npar      = length(lavaan::coef(fit)),
        converged = lavaan::lavInspect(fit, "converged"),
        partable  = lavaan::parameterTable(fit)
      )
    }, error = function(e) {
      NULL
    })
  })

  ok <- vapply(fits, function(f) !is.null(f) && f$converged, logical(1))
  M_ok <- sum(ok)

  if (M_ok < 2) {
    return(list(converged_count = M_ok, success = FALSE))
  }

  fits_ok <- fits[ok]
  npar <- fits_ok[[1]]$npar
  pnames <- names(fits_ok[[1]]$coef)

  coef_mat   <- matrix(NA, nrow = M_ok, ncol = npar)
  vcov_list  <- vector("list", M_ok)
  loglik_vec <- numeric(M_ok)

  for (i in seq_len(M_ok)) {
    coef_mat[i, ] <- fits_ok[[i]]$coef
    vcov_list[[i]] <- fits_ok[[i]]$vcov
    loglik_vec[i]  <- fits_ok[[i]]$loglik
  }
  colnames(coef_mat) <- pnames

  list(
    coefs             = coef_mat,
    vcovs             = vcov_list,
    logliks           = loglik_vec,
    npar              = npar,
    param_names       = pnames,
    converged_count   = M_ok,
    converged_flags   = ok,
    partable_template = fits_ok[[1]]$partable,
    success           = TRUE
  )
}


#' Pool MI Results via Rubin's Rules and Compute W, B, tr(RIV)
#'
#' @param mi_result Output of \code{fit_single_model_mi}.
#' @param config Configuration list (uses \code{ridge_factor}).
#' @return List with \code{theta_bar}, \code{W}, \code{B}, \code{RIV},
#'   \code{tr_RIV}, \code{mean_loglik}, \code{M_ok}.
pool_mi <- function(mi_result, config) {
  M_ok <- mi_result$converged_count
  Q    <- mi_result$npar

  theta_bar <- colMeans(mi_result$coefs)

  W <- Reduce("+", mi_result$vcovs) / M_ok

  centered <- sweep(mi_result$coefs, 2, theta_bar)
  B <- crossprod(centered) / (M_ok - 1)

  W_inv <- tryCatch({
    ridge <- config$ridge_factor * mean(diag(W))
    solve(W + diag(ridge, nrow = Q))
  }, error = function(e) {
    warning("W matrix singular; using pseudoinverse")
    MASS::ginv(W)
  })

  RIV <- (1 + 1 / M_ok) * W_inv %*% B
  tr_RIV <- sum(diag(RIV))

  # Eigenvalue spectrum of RIV via Cholesky-symmetrized eigendecomposition
  # of (1 + 1/M) W^{-1/2}^T B W^{-1/2}. Numerically stable via
  # eigen(symmetric = TRUE) on a symmetric matrix similar to W^{-1} B.
  # Used by v4.5 §13 finite-M variance correction (chi2_MI_corr).
  spec <- tryCatch({
    R_chol <- chol(W + diag(config$ridge_factor * mean(diag(W)), nrow = Q))
    R_inv  <- backsolve(R_chol, diag(Q))
    A_sym  <- crossprod(R_inv, B) %*% R_inv
    A_sym  <- (A_sym + t(A_sym)) / 2  # numerical symmetrization
    ev     <- eigen(A_sym, symmetric = TRUE, only.values = TRUE)
    lambdas <- (1 + 1 / M_ok) * ev$values
    list(lambdas = lambdas, sum_lambda_sq = sum(lambdas^2))
  }, error = function(e) {
    list(lambdas = rep(NA_real_, Q), sum_lambda_sq = NA_real_)
  })

  mean_loglik <- mean(mi_result$logliks)

  list(
    theta_bar     = theta_bar,
    W             = W,
    B             = B,
    RIV           = RIV,
    tr_RIV        = tr_RIV,
    lambdas       = spec$lambdas,
    sum_lambda_sq = spec$sum_lambda_sq,
    mean_loglik   = mean_loglik,
    M_ok          = M_ok
  )
}


#' Build a Parameter Table with Free Parameters Fixed to Given Values
#'
#' @param partable_template Parameter table from a converged fit.
#' @param theta Named numeric vector of parameter values.
#' @return A parameter table with free parameters fixed to theta, or NULL
#'   if any parameter name is missing from theta.
build_fixed_partable <- function(partable_template, theta) {
  partable_fixed <- partable_template
  free_idx <- partable_fixed$free > 0

  param_names <- ifelse(
    partable_fixed$label[free_idx] != "",
    partable_fixed$label[free_idx],
    paste0(partable_fixed$lhs[free_idx],
           partable_fixed$op[free_idx],
           partable_fixed$rhs[free_idx])
  )
  fixed_values <- theta[param_names]
  if (any(is.na(fixed_values))) return(NULL)

  partable_fixed$ustart[free_idx] <- fixed_values
  partable_fixed$lower[free_idx]  <- fixed_values
  partable_fixed$upper[free_idx]  <- fixed_values
  partable_fixed
}


#' Evaluate Log-Likelihood on a Single Dataset at a Fixed Theta
#'
#' @param data A single data.frame.
#' @param partable_fixed A parameter table with all free parameters fixed
#'   (built via \code{\link{build_fixed_partable}}).
#' @return Numeric log-likelihood, NA on failure.
eval_loglik_at_theta <- function(data, partable_fixed) {
  if (is.null(partable_fixed)) return(NA_real_)
  tryCatch({
    fit_fixed <- lavaan::sem(
      model         = partable_fixed,
      data          = data,
      meanstructure = TRUE,
      do.fit        = TRUE,
      se            = "none",
      test          = "standard",
      warn          = FALSE
    )
    as.numeric(lavaan::logLik(fit_fixed))
  }, error = function(e) NA_real_)
}


#' Evaluate Log-Likelihoods on Each Imputed Dataset at Pooled Theta
#'
#' Uses the fixed-parameter technique: fix all free parameters to pooled
#' values via lower=upper bounds, then fit to get the log-likelihood.
#'
#' @param imputed_list List of M data.frames.
#' @param model_syntax lavaan model syntax (character). Currently unused —
#'   kept for signature compatibility.
#' @param theta_bar Named vector of pooled parameter estimates.
#' @param partable_template Parameter table from a converged fit.
#' @param converged_flags Logical vector indicating which imputations converged.
#' @return Numeric vector of length M (NA for failed evaluations).
eval_loglik_at_pooled <- function(imputed_list, model_syntax, theta_bar,
                                  partable_template, converged_flags) {
  M <- length(imputed_list)
  partable_fixed <- build_fixed_partable(partable_template, theta_bar)
  if (is.null(partable_fixed)) {
    warning("Could not match all parameters to pooled estimates")
    return(rep(NA_real_, M))
  }

  vapply(seq_len(M), function(m) {
    if (!converged_flags[m]) return(NA_real_)
    eval_loglik_at_theta(imputed_list[[m]], partable_fixed)
  }, numeric(1))
}


#' Three-Term Bias-Decomposition Log-Likelihoods (Per Model)
#'
#' Evaluates, for each model, the complete-data log-likelihood at two
#' alternative parameter points besides the oracle \eqn{\hat\theta_{com}}:
#' \itemize{
#'   \item \code{loglik_com_at_obs}: \eqn{\ell_{com}(\hat\theta_{obs})}
#'         — complete-data log-lik at the FIML observed-data MLE.
#'   \item \code{loglik_com_at_pooled}: \eqn{\ell_{com}(\bar\theta_{MI})}
#'         — complete-data log-lik at Rubin's-rule pooled estimate.
#' }
#' Combined with the oracle \eqn{\ell_{com}(\hat\theta_{com})} (already in
#' \code{complete_fits[[m]]$loglik}) and the pooled
#' \eqn{\bar Q_{MI}(\bar\theta_{MI})} (mean of
#' \code{mi_fits[[m]]$logliks_at_pooled}), these allow per-rep
#' decomposition of the pooled-MI vs oracle bias into Imputation,
#' Pooling-approximation, and Estimation-mismatch pieces.
#'
#' @param complete_fits Named list from \code{\link{fit_complete}}; needs
#'   \code{partable_template} (stored by \code{fit_complete}).
#' @param observed_fits Named list from \code{\link{fit_observed}}; needs
#'   \code{coef}.
#' @param mi_fits Named list from \code{\link{fit_mi_models}}; needs
#'   \code{theta_bar}.
#' @param data_complete Complete-data data.frame.
#' @return data.frame with columns \code{model},
#'   \code{loglik_com_at_obs}, \code{loglik_com_at_pooled}.
compute_decomp_logliks <- function(complete_fits, observed_fits,
                                   mi_fits, data_complete) {
  model_names <- names(complete_fits)

  rows <- lapply(model_names, function(mname) {
    cf <- complete_fits[[mname]]
    of <- observed_fits[[mname]]
    mf <- mi_fits[[mname]]

    ll_at_obs    <- NA_real_
    ll_at_pooled <- NA_real_

    if (isTRUE(cf$converged) && !is.null(cf$partable_template)) {

      # ell_com(theta_obs_fiml)
      if (isTRUE(of$converged) && !is.null(of$coef)) {
        ptab_obs  <- build_fixed_partable(cf$partable_template, of$coef)
        ll_at_obs <- eval_loglik_at_theta(data_complete, ptab_obs)
      }

      # ell_com(theta_pooled)
      if (isTRUE(mf$success) && !is.null(mf$theta_bar)) {
        ptab_pool    <- build_fixed_partable(cf$partable_template, mf$theta_bar)
        ll_at_pooled <- eval_loglik_at_theta(data_complete, ptab_pool)
      }
    }

    data.frame(
      model                = mname,
      loglik_com_at_obs    = ll_at_obs,
      loglik_com_at_pooled = ll_at_pooled,
      stringsAsFactors     = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Fit All Candidate + Saturated Models to MI Data and Pool
#'
#' For each model: fit to all imputations, pool via Rubin's rules,
#' and evaluate log-likelihoods at pooled estimates.
#'
#' @param imputed_list List of M data.frames.
#' @param models Named list of lavaan model syntaxes.
#' @param config Configuration list.
#' @param pop_starts Optional named list of warm-start coefficient
#'   vectors (see \code{\link{compute_pop_starts}}).
#' @return Named list (one per model) with pooling results +
#'   \code{logliks_at_pooled}.
fit_mi_models <- function(imputed_list, models, config, pop_starts = NULL) {
  results <- lapply(names(models), function(mname) {
    start_coefs <- if (!is.null(pop_starts)) pop_starts[[mname]] else NULL
    mi_raw <- fit_single_model_mi(imputed_list, models[[mname]], start_coefs = start_coefs)

    if (!mi_raw$success) {
      return(list(success = FALSE, model = mname,
                  converged_count = mi_raw$converged_count))
    }

    pooled <- pool_mi(mi_raw, config)

    logliks_at_pooled <- eval_loglik_at_pooled(
      imputed_list      = imputed_list,
      model_syntax      = models[[mname]],
      theta_bar         = pooled$theta_bar,
      partable_template = mi_raw$partable_template,
      converged_flags   = mi_raw$converged_flags
    )

    list(
      success           = TRUE,
      model             = mname,
      npar              = mi_raw$npar,
      theta_bar         = pooled$theta_bar,
      W                 = pooled$W,
      B                 = pooled$B,
      tr_RIV            = pooled$tr_RIV,
      lambdas           = pooled$lambdas,
      sum_lambda_sq     = pooled$sum_lambda_sq,
      mean_loglik       = pooled$mean_loglik,
      logliks           = mi_raw$logliks,
      logliks_at_pooled = logliks_at_pooled,
      converged_count   = mi_raw$converged_count
    )
  })
  names(results) <- names(models)
  return(results)
}
