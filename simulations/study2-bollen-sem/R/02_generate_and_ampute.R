# ============================================================================
# Data Generation and Amputation for Study 2 (SEM)
# ============================================================================

#' Generate Complete Data from Population Covariance
#'
#' @param n Sample size
#' @param config Configuration list from get_config()
#' @param seed Random seed
#' @return data.frame with y1-y9
generate_complete_data <- function(n, config, seed) {
  set.seed(seed)
  dat <- MASS::mvrnorm(n = n, mu = rep(0, 9), Sigma = config$sigma_pop)
  dat <- as.data.frame(dat)
  colnames(dat) <- config$var_names
  return(dat)
}


#' Ampute Data (Introduce MCAR Missingness)
#'
#' @param data Complete data.frame
#' @param miss_rate Proportion of incomplete cases (0-1)
#' @param seed Random seed
#' @return data.frame with NAs
ampute_data <- function(data, miss_rate, seed) {
  set.seed(seed)

  # mice::ampute returns a mads object; $amp is the amputed data
  amp_result <- mice::ampute(
    data    = data,
    prop    = miss_rate,
    mech    = "MCAR"
  )

  return(amp_result$amp)
}
