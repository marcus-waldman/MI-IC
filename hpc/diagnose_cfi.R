suppressPackageStartupMessages({
  library(miicsem); library(lavaan); library(mice)
})
config <- get_config(n_reps = 5, master_seed = 32897891,
                     sample_sizes = 500, miss_rates = 0.4, M = 5)
seeds <- generate_seeds(config$master_seed, 5)
var_names <- config$var_names
m1_syntax  <- get_sim1_models()$M1
sat_syntax <- get_saturated_model(var_names)
null_syntax <- paste(sprintf("%s ~~ %s", var_names, var_names), collapse = "\n")

dc <- miicsem:::generate_complete_data(500, config, seeds[[1]]$data)
cat("dc dim:", dim(dc), "\n")

fit1 <- tryCatch(sem(m1_syntax, data = dc, meanstructure = TRUE, warn = FALSE),
                 error = function(e) { cat("M1 err:", conditionMessage(e), "\n"); NULL })
cat("M1 npar:", length(coef(fit1)), "  loglik:", as.numeric(logLik(fit1)), "\n")

fitsat <- tryCatch(sem(sat_syntax, data = dc, meanstructure = TRUE, warn = FALSE),
                   error = function(e) { cat("Msat err:", conditionMessage(e), "\n"); NULL })
cat("Msat npar:", length(coef(fitsat)), "  loglik:", as.numeric(logLik(fitsat)), "\n")

fit0 <- tryCatch(sem(null_syntax, data = dc, meanstructure = TRUE, warn = FALSE),
                 error = function(e) { cat("M0 err:", conditionMessage(e), "\n"); NULL })
cat("M0 npar:", length(coef(fit0)), "  loglik:", as.numeric(logLik(fit0)), "\n")

theta_M1 <- coef(fit1)
cat("theta_M1 length:", length(theta_M1), "\n")
fit1_at <- tryCatch(sem(m1_syntax, data = dc, meanstructure = TRUE,
                        start = theta_M1, do.fit = FALSE, warn = FALSE),
                    error = function(e) { cat("at err:", conditionMessage(e), "\n"); NULL })

# Try implied moments
implied <- tryCatch(lavInspect(fit1_at, "implied"),
                    error = function(e) { cat("implied err:", conditionMessage(e), "\n"); NULL })
if (!is.null(implied)) {
  cat("implied names:", paste(names(implied), collapse=","), "\n")
  cat("mean[1:3]:", head(implied$mean, 3), "\n")
  cat("cov dim:", dim(implied$cov), "\n")
}

# Compute Gaussian loglik manually
gauss_ll <- function(Y, mu, Sigma) {
  Y <- as.matrix(Y); N <- nrow(Y); q <- ncol(Y)
  ybar <- colMeans(Y)
  S    <- crossprod(scale(Y, scale = FALSE)) / N   # MLE
  invS <- solve(Sigma)
  -0.5 * N * (q * log(2*pi) + log(det(Sigma)) +
              sum(diag(S %*% invS)) +
              as.numeric(t(ybar - mu) %*% invS %*% (ybar - mu)))
}

ll_at <- gauss_ll(dc, implied$mean, implied$cov)
cat("manual loglik at theta_M1 (own MLE):", ll_at, "\n")
cat("lavaan logLik(fit1):                ", as.numeric(logLik(fit1)), "\n")
cat("difference:", ll_at - as.numeric(logLik(fit1)), "\n")

# Quick perturbation test: shift theta and see if loglik decreases
theta_pert <- theta_M1 + 0.01
fit1_pert <- sem(m1_syntax, data = dc, meanstructure = TRUE,
                 start = theta_pert, do.fit = FALSE, warn = FALSE)
imp_pert <- lavInspect(fit1_pert, "implied")
ll_pert <- gauss_ll(dc, imp_pert$mean, imp_pert$cov)
cat("loglik at perturbed theta:           ", ll_pert, "\n")
cat("diff (should be negative):", ll_pert - ll_at, "\n")
