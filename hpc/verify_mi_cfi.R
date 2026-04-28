# ============================================================================
# Standalone verification of the MI-CFI claims in v4.3 Section 10
# (mi_deviance_bias_derivation_v4.qmd).
#
# Predictions to verify (Bollen 3-factor CFA, q=9, p_M1=32, p_M0=18,
# p_Msat=54, df_M1=22, df_M0=36):
#   1. CFI_com - CFI_adhoc ≈ Δ(M1) / [chi2_com(M0) - df_0]
#      where Δ(M1) = 2 · [tr(RIV)_Msat − tr(RIV)_M1].
#      With Δ(M1) ≈ 5 and chi2_com(M0) ≈ 1500–2000, the bias should be
#      ~0.003 (a few thousandths).
#   2. CFI_MI ≈ CFI_com to within MC error.
#
# Design: PMM, N=500, mr ∈ {0.10, 0.25, 0.40}, M=50, 200 reps, 16 cores.
# Runs locally; ~30–60 min total.
#
# Output: hpc/results-mi-cfi/cfi_summary.rds and stdout table.
# ============================================================================

suppressPackageStartupMessages({
  library(miicsem); library(lavaan); library(mice)
  library(parallel); library(pbapply)
})

cat(sprintf("miicsem version: %s\n", as.character(packageVersion("miicsem"))))

N           <- 500L
MR_VALUES   <- c(0.10, 0.25, 0.40)
M_IMP       <- 50L
N_REPS      <- 200L
N_CORES     <- 16L
BASE_SEED   <- 32897891L

results_dir <- "hpc/results-mi-cfi"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

config <- get_config(n_reps = N_REPS, master_seed = BASE_SEED,
                     sample_sizes = N, miss_rates = MR_VALUES, M = M_IMP)
seeds <- generate_seeds(config$master_seed, config$n_reps)
var_names <- config$var_names

m1_syntax  <- get_sim1_models()$M1
sat_syntax <- get_saturated_model(var_names)
null_syntax <- paste(sprintf("%s ~~ %s", var_names, var_names), collapse = "\n")

# --- Helpers ---------------------------------------------------------------
# Standard SEM fit; returns NULL on failure or non-convergence.
fit_lav <- function(syntax, data) {
  tryCatch({
    fit <- sem(syntax, data = data, meanstructure = TRUE, warn = FALSE)
    if (!lavInspect(fit, "converged")) return(NULL)
    list(loglik = as.numeric(logLik(fit)),
         coef   = coef(fit),
         vcov   = vcov(fit),
         npar   = length(coef(fit)))
  }, error = function(e) NULL)
}

# Gaussian log-likelihood at a given (mu, Sigma) for complete data.
gauss_ll <- function(Y, mu, Sigma) {
  Y <- as.matrix(Y); N <- nrow(Y); q <- ncol(Y)
  ybar <- colMeans(Y)
  S    <- crossprod(scale(Y, scale = FALSE)) / N
  invS <- solve(Sigma)
  -0.5 * N * (q * log(2*pi) + log(det(Sigma)) +
              sum(diag(S %*% invS)) +
              as.numeric(t(ybar - mu) %*% invS %*% (ybar - mu)))
}

# Loglik of `data` at the model implied by `theta` for `syntax`. Uses lavaan
# do.fit=FALSE to read off implied (mu, Sigma), then evaluates Gaussian ll
# manually (lavaan's logLik() requires actual estimation, which we can't do
# at a fixed point).
ll_at_theta <- function(syntax, data, theta) {
  tryCatch({
    fit <- sem(syntax, data = data, meanstructure = TRUE,
               start = theta, do.fit = FALSE, warn = FALSE)
    imp <- lavInspect(fit, "implied")
    gauss_ll(data, imp$mean, imp$cov)
  }, error = function(e) NA_real_)
}

run_one_cfi_rep <- function(rep_id, n_val, mr_val) {
  s <- seeds[[rep_id]]

  # Complete data + complete fits
  dc <- miicsem:::generate_complete_data(n_val, config, s$data)
  fit_M1_com  <- fit_lav(m1_syntax,   dc)
  fit_sat_com <- fit_lav(sat_syntax,  dc)
  fit_M0_com  <- fit_lav(null_syntax, dc)
  if (is.null(fit_M1_com) || is.null(fit_sat_com) || is.null(fit_M0_com))
    return(NULL)

  DEV_com_M1  <- -2 * fit_M1_com$loglik
  DEV_com_sat <- -2 * fit_sat_com$loglik
  DEV_com_M0  <- -2 * fit_M0_com$loglik
  chi2_com_M1 <- DEV_com_M1 - DEV_com_sat
  chi2_com_M0 <- DEV_com_M0 - DEV_com_sat

  # Ampute + impute
  dm <- miicsem:::ampute_data(dc, mr_val, s$ampute)
  imp <- tryCatch(
    mice::mice(dm, m = config$M, method = "pmm",
               maxit = config$mice_maxit, seed = s$impute, printFlag = FALSE),
    error = function(e) NULL)
  if (is.null(imp)) return(NULL)
  imps <- lapply(seq_len(config$M), function(m) mice::complete(imp, action = m))

  # Per-imputation own-MLE fits
  fits_M1  <- lapply(imps, function(d) fit_lav(m1_syntax,   d))
  fits_sat <- lapply(imps, function(d) fit_lav(sat_syntax,  d))
  fits_M0  <- lapply(imps, function(d) fit_lav(null_syntax, d))
  if (any(vapply(fits_M1,  is.null, TRUE)) ||
      any(vapply(fits_sat, is.null, TRUE)) ||
      any(vapply(fits_M0,  is.null, TRUE))) return(NULL)

  # Pool: bar_theta, tr(RIV), DEV_adhoc per model
  pool_one <- function(fits) {
    coefs <- do.call(cbind, lapply(fits, function(f) f$coef))
    bar_theta <- rowMeans(coefs)
    W <- Reduce("+", lapply(fits, function(f) f$vcov)) / length(fits)
    B <- stats::cov(t(coefs))
    tr_RIV <- (1 + 1/length(fits)) * sum(diag(solve(W) %*% B))
    list(bar_theta = bar_theta, npar = length(bar_theta), tr_RIV = tr_RIV,
         own_ll = vapply(fits, function(f) f$loglik, numeric(1)))
  }
  pooled_M1  <- pool_one(fits_M1)
  pooled_sat <- pool_one(fits_sat)
  pooled_M0  <- pool_one(fits_M0)

  # Q-bar at pooled theta: per-imputation loglik at bar_theta
  ll_M1_pooled  <- vapply(imps, function(d) ll_at_theta(m1_syntax,   d, pooled_M1$bar_theta),  numeric(1))
  ll_sat_pooled <- vapply(imps, function(d) ll_at_theta(sat_syntax,  d, pooled_sat$bar_theta), numeric(1))
  ll_M0_pooled  <- vapply(imps, function(d) ll_at_theta(null_syntax, d, pooled_M0$bar_theta),  numeric(1))
  if (any(is.na(ll_M1_pooled)) || any(is.na(ll_sat_pooled)) ||
      any(is.na(ll_M0_pooled))) return(NULL)

  DEV_adhoc_M1  <- -2 * mean(pooled_M1$own_ll)
  DEV_adhoc_sat <- -2 * mean(pooled_sat$own_ll)
  DEV_adhoc_M0  <- -2 * mean(pooled_M0$own_ll)
  MI_DEV_M1  <- -2 * mean(ll_M1_pooled)  + pooled_M1$tr_RIV
  MI_DEV_sat <- -2 * mean(ll_sat_pooled) + pooled_sat$tr_RIV
  MI_DEV_M0  <- -2 * mean(ll_M0_pooled)  + pooled_M0$tr_RIV

  chi2_adhoc_M1 <- DEV_adhoc_M1 - DEV_adhoc_sat
  chi2_adhoc_M0 <- DEV_adhoc_M0 - DEV_adhoc_sat
  chi2_MI_M1    <- MI_DEV_M1    - MI_DEV_sat
  chi2_MI_M0    <- MI_DEV_M0    - MI_DEV_sat

  df_M1 <- pooled_sat$npar - pooled_M1$npar
  df_M0 <- pooled_sat$npar - pooled_M0$npar

  cfi_one <- function(chi_M, df_M, chi_0, df_0) {
    num <- max(chi_M - df_M, 0); den <- max(chi_M - df_M, chi_0 - df_0, 0)
    if (den <= 0) return(1.0)
    1 - num / den
  }
  cfi_com   <- cfi_one(chi2_com_M1,   df_M1, chi2_com_M0,   df_M0)
  cfi_adhoc <- cfi_one(chi2_adhoc_M1, df_M1, chi2_adhoc_M0, df_M0)
  cfi_MI    <- cfi_one(chi2_MI_M1,    df_M1, chi2_MI_M0,    df_M0)

  data.frame(
    rep_id = rep_id, n = n_val, mr = mr_val,
    df_M1 = df_M1, df_M0 = df_M0,
    chi2_com_M1   = chi2_com_M1,   chi2_com_M0   = chi2_com_M0,
    chi2_adhoc_M1 = chi2_adhoc_M1, chi2_adhoc_M0 = chi2_adhoc_M0,
    chi2_MI_M1    = chi2_MI_M1,    chi2_MI_M0    = chi2_MI_M0,
    tr_RIV_M1 = pooled_M1$tr_RIV, tr_RIV_M0 = pooled_M0$tr_RIV,
    tr_RIV_sat = pooled_sat$tr_RIV,
    cfi_com = cfi_com, cfi_adhoc = cfi_adhoc, cfi_MI = cfi_MI
  )
}

# --- Driver -----------------------------------------------------------------
overall_t0 <- proc.time()

all_summary <- list()
for (mr in MR_VALUES) {
  cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
  out_file <- file.path(results_dir, sprintf("cfi_%s.rds", cell_label))
  cat(sprintf("\n=== %s ===\n", cell_label))
  if (file.exists(out_file)) {
    cat("  cached, loading\n")
    rep_rows <- readRDS(out_file)
  } else {
    t0 <- proc.time()
    cl <- makeCluster(N_CORES)
    clusterExport(cl,
      varlist = c("seeds", "config", "var_names",
                  "m1_syntax", "sat_syntax", "null_syntax",
                  "fit_lav", "ll_at_theta", "gauss_ll", "run_one_cfi_rep",
                  "N", "mr"),
      envir = environment())
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(miicsem); library(lavaan); library(mice)
      })
    })
    rep_rows <- pblapply(seq_len(N_REPS), function(i) {
      tryCatch(run_one_cfi_rep(i, N, mr), error = function(e) NULL)
    }, cl = cl)
    stopCluster(cl)
    n_ok <- sum(!vapply(rep_rows, is.null, TRUE))
    rep_rows <- if (n_ok > 0) do.call(rbind, rep_rows[!vapply(rep_rows, is.null, TRUE)]) else NULL
    if (!is.null(rep_rows)) saveRDS(rep_rows, out_file)
    cat(sprintf("  %d/%d reps OK in %.1fs\n",
                n_ok, N_REPS, (proc.time() - t0)[3]))
  }
  all_summary[[cell_label]] <- rep_rows
}

saveRDS(all_summary, file.path(results_dir, "cfi_summary.rds"))
cat(sprintf("\n\nWall time: %.1fs (%.2fh)\n",
            (proc.time() - overall_t0)[3],
            (proc.time() - overall_t0)[3] / 3600))

cat("\n=== Empirical CFI variants vs analytic prediction ===\n")
cat(sprintf("%-6s %5s | %8s %10s %8s | %8s %10s %8s | %10s %12s %8s\n",
            "mr", "n_rep",
            "chi2com1", "chi2adhoc1", "chi2MI1",
            "chi2com0", "chi2adhoc0", "chi2MI0",
            "CFI_com", "CFI_adhoc", "CFI_MI"))
for (mr in MR_VALUES) {
  cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
  d <- all_summary[[cell_label]]
  if (is.null(d) || nrow(d) == 0) next
  cat(sprintf("%-6.2f %5d | %8.2f %10.2f %8.2f | %8.1f %10.1f %8.1f | %10.4f %12.4f %8.4f\n",
              mr, nrow(d),
              mean(d$chi2_com_M1), mean(d$chi2_adhoc_M1), mean(d$chi2_MI_M1),
              mean(d$chi2_com_M0), mean(d$chi2_adhoc_M0), mean(d$chi2_MI_M0),
              mean(d$cfi_com),     mean(d$cfi_adhoc),     mean(d$cfi_MI)))
}

cat("\nBias check (CFI_com − CFI_x) and analytic prediction:\n")
cat(sprintf("%-6s | %12s %12s | %12s %12s\n",
            "mr", "bias_adhoc", "bias_MI", "predicted", "ratio"))
for (mr in MR_VALUES) {
  cell_label <- sprintf("n=%d_mr=%.2f", N, mr)
  d <- all_summary[[cell_label]]
  if (is.null(d) || nrow(d) == 0) next
  bias_adhoc <- mean(d$cfi_com) - mean(d$cfi_adhoc)
  bias_MI    <- mean(d$cfi_com) - mean(d$cfi_MI)
  delta_M1   <- 2 * (mean(d$tr_RIV_sat) - mean(d$tr_RIV_M1))
  null_disc  <- mean(d$chi2_com_M0) - mean(d$df_M0)
  predicted  <- delta_M1 / null_disc
  ratio      <- if (abs(predicted) > 1e-6) bias_adhoc / predicted else NA
  cat(sprintf("%-6.2f | %12.5f %12.5f | %12.5f %12.3f\n",
              mr, bias_adhoc, bias_MI, predicted, ratio))
}

cat("\nDone.\n")
