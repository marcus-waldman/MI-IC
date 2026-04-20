#' Fit All 12 Models to Complete Data (Oracle)
#'
#' @param data_complete Complete data.frame (no missingness).
#' @param models Named list of lavaan model syntaxes.
#' @return Named list; each element has \code{loglik}, \code{npar},
#'   \code{converged}.
fit_complete <- function(data_complete, models) {
  results <- lapply(names(models), function(mname) {
    tryCatch({
      fit <- lavaan::sem(
        model     = models[[mname]],
        data      = data_complete,
        estimator = "ML",
        se        = "none",
        test      = "standard",
        warn      = FALSE
      )
      list(
        loglik    = as.numeric(lavaan::logLik(fit)),
        npar      = lavaan::fitmeasures(fit, "npar"),
        converged = lavaan::lavInspect(fit, "converged")
      )
    }, error = function(e) {
      list(loglik = NA, npar = NA, converged = FALSE)
    })
  })
  names(results) <- names(models)
  return(results)
}


#' Fit a Single Model to All M Imputed Datasets
#'
#' @param imputed_list List of M data.frames.
#' @param model_syntax lavaan model syntax (character).
#' @return List with \code{coefs}, \code{vcovs}, \code{logliks},
#'   \code{npar}, \code{param_names}, \code{converged_count},
#'   \code{converged_flags}, \code{partable_template}, \code{success}.
fit_single_model_mi <- function(imputed_list, model_syntax) {
  M <- length(imputed_list)

  fits <- lapply(seq_len(M), function(m) {
    tryCatch({
      fit <- lavaan::sem(
        model     = model_syntax,
        data      = imputed_list[[m]],
        estimator = "ML",
        se        = "standard",
        test      = "standard",
        warn      = FALSE
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

  mean_loglik <- mean(mi_result$logliks)

  list(
    theta_bar   = theta_bar,
    W           = W,
    B           = B,
    RIV         = RIV,
    tr_RIV      = tr_RIV,
    mean_loglik = mean_loglik,
    M_ok        = M_ok
  )
}


#' Evaluate Log-Likelihoods on Each Imputed Dataset at Pooled Theta
#'
#' Uses the fixed-parameter technique: fix all free parameters to pooled
#' values via lower=upper bounds, then fit to get the log-likelihood.
#'
#' @param imputed_list List of M data.frames.
#' @param model_syntax lavaan model syntax (character).
#' @param theta_bar Named vector of pooled parameter estimates.
#' @param partable_template Parameter table from a converged fit.
#' @param converged_flags Logical vector indicating which imputations converged.
#' @return Numeric vector of length M (NA for failed evaluations).
eval_loglik_at_pooled <- function(imputed_list, model_syntax, theta_bar,
                                  partable_template, converged_flags) {
  M <- length(imputed_list)

  partable_fixed <- partable_template
  free_idx <- partable_fixed$free > 0

  param_names <- ifelse(
    partable_fixed$label[free_idx] != "",
    partable_fixed$label[free_idx],
    paste0(partable_fixed$lhs[free_idx],
           partable_fixed$op[free_idx],
           partable_fixed$rhs[free_idx])
  )
  pooled_values <- theta_bar[param_names]

  if (any(is.na(pooled_values))) {
    warning("Could not match all parameters to pooled estimates")
    return(rep(NA, M))
  }

  partable_fixed$ustart[free_idx] <- pooled_values
  partable_fixed$lower[free_idx]  <- pooled_values
  partable_fixed$upper[free_idx]  <- pooled_values

  logliks <- vapply(seq_len(M), function(m) {
    if (!converged_flags[m]) return(NA_real_)
    tryCatch({
      fit_fixed <- lavaan::sem(
        model  = partable_fixed,
        data   = imputed_list[[m]],
        do.fit = TRUE,
        se     = "none",
        test   = "standard",
        warn   = FALSE
      )
      as.numeric(lavaan::logLik(fit_fixed))
    }, error = function(e) {
      NA_real_
    })
  }, numeric(1))

  return(logliks)
}


#' Fit All 12 Models to MI Data and Pool
#'
#' For each model: fit to all imputations, pool via Rubin's rules,
#' and evaluate log-likelihoods at pooled estimates.
#'
#' @param imputed_list List of M data.frames.
#' @param models Named list of lavaan model syntaxes.
#' @param config Configuration list.
#' @return Named list (one per model) with pooling results +
#'   \code{logliks_at_pooled}.
fit_mi_models <- function(imputed_list, models, config) {
  results <- lapply(names(models), function(mname) {
    mi_raw <- fit_single_model_mi(imputed_list, models[[mname]])

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
      mean_loglik       = pooled$mean_loglik,
      logliks_at_pooled = logliks_at_pooled,
      converged_count   = mi_raw$converged_count
    )
  })
  names(results) <- names(models)
  return(results)
}
