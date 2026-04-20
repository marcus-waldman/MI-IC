#' Generate Complete Data from the Population Covariance
#'
#' @param n Sample size.
#' @param config Configuration list from \code{\link{get_config}}.
#' @param seed Random seed.
#' @return A data.frame with columns y1 through y9.
generate_complete_data <- function(n, config, seed) {
  set.seed(seed)
  dat <- MASS::mvrnorm(n = n, mu = rep(0, 9), Sigma = config$sigma_pop)
  dat <- as.data.frame(dat)
  colnames(dat) <- config$var_names
  return(dat)
}


#' Introduce MCAR Missingness via mice::ampute
#'
#' @param data Complete data.frame.
#' @param miss_rate Proportion of incomplete cases (0-1).
#' @param seed Random seed.
#' @return A data.frame with NAs.
ampute_data <- function(data, miss_rate, seed) {
  set.seed(seed)
  amp_result <- mice::ampute(
    data = data,
    prop = miss_rate,
    mech = "MCAR"
  )
  return(amp_result$amp)
}
