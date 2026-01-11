# ============================================================================
# Data Generation for Study 1
# ============================================================================

#' Generate Mediation Data
#'
#' @param n Sample size
#' @param a Effect of X on M (default: 0.5)
#' @param b Effect of M on Y (default: 0.5)
#' @param sigma2_M Residual variance for M (default: 1)
#' @param sigma2_Y Residual variance for Y (default: 1)
#' @param seed Random seed (optional)
#' @return data.frame with X, M, Y
generate_data <- function(n, a = 0.5, b = 0.5, sigma2_M = 1, sigma2_Y = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- rnorm(n)
  M <- a * X + rnorm(n, 0, sqrt(sigma2_M))
  Y <- b * M + rnorm(n, 0, sqrt(sigma2_Y))

  data.frame(X = X, M = M, Y = Y)
}
