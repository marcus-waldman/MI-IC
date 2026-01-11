# ============================================================================
# Configuration for Study 1: Pedagogical Simulation
# ============================================================================
# All adjustable parameters for the simulation
# ============================================================================

config <- list(

  # --------------------------------------------------------------------------
  # Sample and Missingness
  # --------------------------------------------------------------------------
  n = 100,                    # Sample size

  prop_missing = 0.5,         # Proportion of observations with missingness
  n_reps = 1000,               # Number of replications

  # --------------------------------------------------------------------------
  # Parallel Processing
  # --------------------------------------------------------------------------
  n_cores = 25,             # Number of cores (NULL = detectCores() - 1)
  parallel_type = "PSOCK",    # Cluster type: "PSOCK" (all OS) or "FORK" (Unix only)

  # --------------------------------------------------------------------------
  # Random Seeds
  # --------------------------------------------------------------------------
  master_seed = 20250110,     # Master seed for SeedMaker
  seeds_file = "results/seeds.rds",  # Path to save/load pre-generated seeds

  # --------------------------------------------------------------------------
  # brms Settings
  # --------------------------------------------------------------------------
  brms_backend = "rstan",     # "rstan" or "cmdstanr"
  brms_chains = 1,            # Number of chains (1 for parallelization across reps)
  brms_warmup = 1000,         # Warmup iterations
  brms_iter = 3000,           # Total iterations (post-warmup = iter - warmup = 2000)
  brms_refresh = 0,           # Progress output (0 = silent)

  # --------------------------------------------------------------------------
  # Population Parameters (DGM: Full Mediation)
  # --------------------------------------------------------------------------
  # X -> M -> Y with no direct X -> Y effect
  beta_xm = 0.5,              # X -> M effect

  beta_my = 0.5,              # M -> Y effect
  beta_xy = 0,                # X -> Y direct effect (0 = full mediation)
  sigma2_m = 1,               # Residual variance for M
  sigma2_y = 1,               # Residual variance for Y
  var_x = 1,                  # Variance of X

  # --------------------------------------------------------------------------
  # Factorial Design
  # --------------------------------------------------------------------------
  # These define the 3x3 conditions
  mechanisms = c("MCAR", "MAR", "MNAR"),
  imputation_models = c("matched", "saturated", "uncongenial"),

  # --------------------------------------------------------------------------
  # Seed Offsets for Conditions
  # --------------------------------------------------------------------------
  # Within a replication, different conditions use seed + offset
  # This ensures reproducibility while allowing condition-specific randomness
  mechanism_offsets = c(MCAR = 0, MAR = 1000, MNAR = 2000),
  model_offsets = c(matched = 0, saturated = 100, uncongenial = 200)
)


# ============================================================================
# Derived Settings (computed from config)
# ============================================================================

#' Get Full Condition Grid
#'
#' @return data.frame with all 9 conditions
get_condition_grid <- function(config) {
  expand.grid(
    mechanism = config$mechanisms,
    imputation_model = config$imputation_models,
    stringsAsFactors = FALSE
  )
}


#' Get Number of Cores to Use
#'
#' @return Integer number of cores
get_n_cores <- function(config) {
  if (is.null(config$n_cores)) {
    max(1, parallel::detectCores() - 1)
  } else {
    config$n_cores
  }
}


#' Get Population Parameters as List
#'
#' @return List with all population parameters
get_pop_params <- function(config) {
  list(
    beta_xm = config$beta_xm,
    beta_my = config$beta_my,
    beta_xy = config$beta_xy,
    sigma2_m = config$sigma2_m,
    sigma2_y = config$sigma2_y,
    var_x = config$var_x,
    indirect_effect = config$beta_xm * config$beta_my,
    direct_effect = config$beta_xy,
    total_effect = config$beta_xm * config$beta_my + config$beta_xy
  )
}
