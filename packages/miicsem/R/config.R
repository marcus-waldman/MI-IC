#' Get Default Simulation Configuration
#'
#' Returns the default design grid and population covariance matrix
#' for the Bollen et al. (2014) SIM1 simulation.
#'
#' @param n_reps Number of replications per condition (default 1000).
#' @param master_seed Master seed for reproducibility (default 32897891).
#' @param sample_sizes Numeric vector of sample sizes.
#' @param miss_rates Numeric vector of missing rates.
#' @param M Number of imputations (default 20).
#'
#' @return A list with all configuration elements.
#' @export
get_config <- function(n_reps = 1000L,
                       master_seed = 32897891L,
                       sample_sizes = c(100, 250, 500, 1000, 5000),
                       miss_rates = c(0.10, 0.25, 0.40),
                       M = 20L) {
  list(
    sample_sizes = sample_sizes,
    miss_rates   = miss_rates,
    miss_mech    = "MCAR",

    M           = as.integer(M),
    mice_method = "pmm",
    mice_maxit  = 10L,

    n_reps      = as.integer(n_reps),
    master_seed = as.integer(master_seed),

    ridge_factor = 1e-6,

    sigma_pop = matrix(c(
      1.0000, 0.4900, 0.4900, 0.4410, 0.2940, 0.3469, 0.2646, 0.1764, 0.1764,
      0.4900, 1.0000, 0.4900, 0.4410, 0.2940, 0.3469, 0.2646, 0.1764, 0.1764,
      0.4900, 0.4900, 1.0000, 0.4410, 0.2940, 0.3469, 0.2646, 0.1764, 0.1764,
      0.4410, 0.4410, 0.4410, 1.0000, 0.5782, 0.6823, 0.5204, 0.3469, 0.3469,
      0.2940, 0.2940, 0.2940, 0.5782, 1.0000, 0.5782, 0.4410, 0.2940, 0.2940,
      0.3469, 0.3469, 0.3469, 0.6823, 0.5782, 1.0000, 0.6145, 0.4410, 0.4410,
      0.2646, 0.2646, 0.2646, 0.5204, 0.4410, 0.6145, 1.0000, 0.5782, 0.5782,
      0.1764, 0.1764, 0.1764, 0.3469, 0.2940, 0.4410, 0.5782, 1.0000, 0.4900,
      0.1764, 0.1764, 0.1764, 0.3469, 0.2940, 0.4410, 0.5782, 0.4900, 1.0000
    ), nrow = 9, ncol = 9, byrow = TRUE,
    dimnames = list(paste0("y", 1:9), paste0("y", 1:9))),

    var_names = paste0("y", 1:9),

    model_categories = list(
      true            = "M1",
      overspecified   = c("M6", "M8", "M9"),
      underspecified  = c("M2", "M3", "M4", "M10"),
      mixed           = c("M5", "M7"),
      wrong_structure = c("M11", "M12")
    ),

    ic_methods = c("AIC_com", "BIC_com",
                   "AIC_adhoc", "BIC_adhoc",
                   "AICcd", "MI_AIC", "MI_BIC")
  )
}
