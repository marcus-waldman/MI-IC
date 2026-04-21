#' Compute Warm-Start Coefficients from the Population Covariance
#'
#' Fits each lavaan model to the population covariance matrix with a very
#' large effective sample size, producing "infinite-N projection" parameter
#' values. The returned named coefficient vectors can be passed via the
#' \code{start} argument to subsequent \code{lavaan::sem()} calls to speed
#' up convergence on imputed datasets.
#'
#' Starting values do not change where lavaan converges (the MLE is the
#' MLE), only how many iterations it takes to get there.
#'
#' @param sigma_pop Population covariance matrix (p x p).
#' @param models Named list of lavaan model syntaxes (as produced by
#'   \code{\link{get_sim1_models}} and \code{\link{get_saturated_model}}).
#' @param nobs_pop Effective sample size used when fitting to
#'   \code{sigma_pop}. Default 1e6 — large enough that the fit is
#'   essentially the exact projection.
#'
#' @return Named list with one element per model; each element is either
#'   a named numeric vector of starting coefficients, or \code{NULL} if
#'   the population fit failed to converge.
compute_pop_starts <- function(sigma_pop, models, nobs_pop = 1e6) {
  pop_starts <- lapply(names(models), function(mname) {
    tryCatch({
      fit <- lavaan::sem(
        model         = models[[mname]],
        sample.cov    = sigma_pop,
        sample.mean   = rep(0, nrow(sigma_pop)),
        sample.nobs   = nobs_pop,
        meanstructure = TRUE,
        estimator     = "ML",
        se            = "none",
        test          = "standard",
        warn          = FALSE
      )
      if (!lavaan::lavInspect(fit, "converged")) return(NULL)
      lavaan::coef(fit)
    }, error = function(e) NULL)
  })
  names(pop_starts) <- names(models)
  return(pop_starts)
}
