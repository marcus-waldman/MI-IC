# ============================================================================
# Step 3a (cozy-imagining-waterfall.md): structural check of v4.5 formula
# ============================================================================
# v4.5 §13 closed-form prediction:
#   Var[chi2_MI] = Var[chi2_com] + 2 * sum_j lambda_j^2 + O(1/M)
# where {lambda_j} are eigenvalues of RIV_perp (per-direction MIFs in the
# orthogonal complement of Theta in Phi).
#
# We don't have per-rep sum_j lambda_j^2 saved -- only tr(RIV)_M and
# tr(RIV)_Msat. tr(RIV)_perp = tr(RIV)_Msat - tr(RIV)_M (additive).
# By Cauchy-Schwarz:
#   (tr(RIV_perp))^2 / Delta_p  <=  sum_j lambda_j^2  <=  (tr(RIV_perp))^2
#
# Define the empirical "concentration ratio":
#   kappa_hat = (Var[chi2_MI] - Var[chi2_com]) / 2 / (tr(RIV)_Msat - tr(RIV)_M)^2
#
# Structural tests:
#   T1. kappa_hat in [1/Delta_p, 1.0]  (Cauchy-Schwarz bounds)
#   T2. kappa_hat approximately constant across N at fixed (mr, imputer, M)
#       -- spectrum is missing-info-pattern-driven, not N-driven.
#   T3. M-scaling: Var[chi2_MI](M=50) - Var[chi2_MI](M=100) ~ B/100 -- O(1/M).
#   T4. Imputer-agnosticism: PMM and Amelia (empri=0) give similar kappa_hat.
#
# Output: diagnostic table + verdict.
# ============================================================================

df <- readRDS("hpc/var_chi2_mi_summary.rds")

cat("=== Var[chi2_MI] aggregate summary loaded ===\n")
cat(sprintf("  rows: %d\n", nrow(df)))
cat(sprintf("  imputers: %s\n", paste(unique(df$imputer), collapse = ", ")))
cat(sprintf("  M_imp: %s\n", paste(sort(unique(df$M_imp)), collapse = ", ")))
cat(sprintf("  N: %s\n", paste(sort(unique(df$n)), collapse = ", ")))
cat(sprintf("  mr: %s\n", paste(sort(unique(df$miss_rate)), collapse = ", ")))
cat(sprintf("  models: %s\n", paste(sort(unique(df$model)), collapse = ", ")))

# Restrict to true model M1 (where mean is calibrated; chi-squares are central)
m1 <- subset(df, model == "M1")

# Implied sum_lambda_sq under v4.5 formula
m1$tr_RIV_perp <- m1$mean_tr_RIV_Msat - m1$mean_tr_RIV_M
m1$Delta_p     <- m1$p_Msat - m1$p_M
m1$emp_excess  <- m1$var_chi2_MI - m1$var_chi2_com
m1$implied_sum_lambda_sq <- m1$emp_excess / 2
m1$kappa_hat   <- m1$implied_sum_lambda_sq / m1$tr_RIV_perp^2
m1$cs_lower    <- 1 / m1$Delta_p   # Cauchy-Schwarz lower bound
m1$cs_upper    <- 1.0
m1$within_CS   <- m1$kappa_hat >= m1$cs_lower & m1$kappa_hat <= m1$cs_upper

# ---- T1: Cauchy-Schwarz bounds ----
cat("\n=== T1: Cauchy-Schwarz check (kappa_hat in [1/Delta_p, 1]) ===\n")
out <- m1[, c("imputer","M_imp","n","miss_rate","n_reps",
              "var_chi2_MI","var_chi2_com","emp_excess",
              "tr_RIV_perp","implied_sum_lambda_sq","kappa_hat",
              "cs_lower","within_CS")]
out <- out[order(out$imputer, out$M_imp, out$n, out$miss_rate), ]
out_round <- out
for (n in names(out_round)) if (is.numeric(out_round[[n]])) out_round[[n]] <- round(out_round[[n]], 3)
print(out_round, row.names = FALSE)

n_pass_CS <- sum(out$within_CS)
cat(sprintf("\n  T1 result: %d/%d cells satisfy CS bounds.\n",
            n_pass_CS, nrow(out)))

# ---- T2: kappa_hat stability across N at fixed (mr, imputer, M) ----
cat("\n=== T2: kappa_hat stability across N at fixed (imputer, M_imp, mr) ===\n")
group_key <- paste(m1$imputer, m1$M_imp, m1$miss_rate, sep = "|")
groups <- split(m1, group_key)

t2_rows <- list()
for (g in names(groups)) {
  s <- groups[[g]]
  s <- s[order(s$n), ]
  if (nrow(s) >= 2) {
    t2_rows[[g]] <- data.frame(
      imputer = s$imputer[1], M_imp = s$M_imp[1], mr = s$miss_rate[1],
      n_cells = nrow(s),
      n_range = paste(range(s$n), collapse = "-"),
      kappa_min = round(min(s$kappa_hat), 3),
      kappa_max = round(max(s$kappa_hat), 3),
      kappa_mean = round(mean(s$kappa_hat), 3),
      kappa_cv = round(sd(s$kappa_hat) / mean(s$kappa_hat), 3),
      stable = abs(max(s$kappa_hat) - min(s$kappa_hat)) < 0.20,  # within 0.2 of each other
      stringsAsFactors = FALSE
    )
  }
}
t2 <- do.call(rbind, t2_rows)
print(t2, row.names = FALSE)
n_stable <- sum(t2$stable, na.rm = TRUE)
cat(sprintf("\n  T2 result: %d/%d (imputer, M, mr) groups have kappa stable across N (range < 0.2).\n",
            n_stable, nrow(t2)))

# ---- T3: M-scaling at fixed (imputer, n, mr) ----
cat("\n=== T3: M-scaling Var[chi2_MI](M=50) vs Var[chi2_MI](M=100) ===\n")
both_M <- subset(m1, imputer == "PMM" &
                     paste(n, miss_rate) %in% intersect(
                       paste(m1$n[m1$M_imp == 50],  m1$miss_rate[m1$M_imp == 50]),
                       paste(m1$n[m1$M_imp == 100], m1$miss_rate[m1$M_imp == 100])))
t3_rows <- list()
for (key in unique(paste(both_M$n, both_M$miss_rate))) {
  s50  <- both_M[both_M$M_imp == 50  & paste(both_M$n, both_M$miss_rate) == key, ]
  s100 <- both_M[both_M$M_imp == 100 & paste(both_M$n, both_M$miss_rate) == key, ]
  if (nrow(s50) == 1 && nrow(s100) == 1) {
    diff <- s50$var_chi2_MI - s100$var_chi2_MI
    # If formula = A + B/M: diff = B/50 - B/100 = B/100, so B = 100*diff.
    # Implied A = Var[chi2_MI](M=100) - B/100 = s100$var - diff
    B_hat <- 100 * diff
    A_hat <- s100$var_chi2_MI - diff
    t3_rows[[key]] <- data.frame(
      n = s50$n, mr = s50$miss_rate,
      var_M50 = s50$var_chi2_MI, var_M100 = s100$var_chi2_MI,
      diff = diff, B_hat = B_hat, A_hat = A_hat,
      var_chi2_com = s100$var_chi2_com,
      excess_at_Minf_implied = A_hat - s100$var_chi2_com,
      stringsAsFactors = FALSE
    )
  }
}
t3 <- do.call(rbind, t3_rows)
t3_round <- t3
for (n in names(t3_round)) if (is.numeric(t3_round[[n]])) t3_round[[n]] <- round(t3_round[[n]], 2)
print(t3_round, row.names = FALSE)
cat(sprintf("\n  T3 interpretation: B_hat is the O(1/M) finite-M coefficient; A_hat = Var[chi2_MI](M=infinity).\n"))
cat(sprintf("  Excess at M=infinity is the v4.5 formula's 2*sum(lambda^2) prediction.\n"))

# ---- T4: imputer-agnosticism ----
cat("\n=== T4: imputer-agnosticism (PMM vs amelia_empri0 at N=500, mr in {0.10,0.25,0.40}, M=50) ===\n")
t4_pmm    <- subset(m1, imputer == "PMM"           & M_imp == 50 & n == 500)
t4_amelia <- subset(m1, imputer == "amelia_empri0" & M_imp == 50 & n == 500)
# Note: PMM at N=500 M=50 may not exist in the grid -- check
if (nrow(t4_pmm) == 0) {
  cat("  No PMM cells at N=500 M=50. Falling back to PMM N=500 M=100 vs amelia N=500 M=50.\n")
  t4_pmm <- subset(m1, imputer == "PMM" & M_imp == 100 & n == 500)
}
cat("\n  PMM:\n")
print(t4_pmm[, c("imputer","M_imp","n","miss_rate","var_chi2_MI","tr_RIV_perp",
                  "implied_sum_lambda_sq","kappa_hat")], row.names = FALSE)
cat("\n  Amelia empri=0:\n")
print(t4_amelia[, c("imputer","M_imp","n","miss_rate","var_chi2_MI","tr_RIV_perp",
                     "implied_sum_lambda_sq","kappa_hat")], row.names = FALSE)

# ---- Save ----
saveRDS(out, "hpc/fit_var_correction_summary.rds")
write.csv(out, "hpc/fit_var_correction_summary.csv", row.names = FALSE)
cat("\n[saved] hpc/fit_var_correction_summary.{rds,csv}\n")

# ---- Verdict ----
cat("\n=== Step 3a verdict ===\n")
verdict_T1 <- if (n_pass_CS == nrow(out)) "PASS" else sprintf("FAIL (%d/%d violate CS)", nrow(out)-n_pass_CS, nrow(out))
cat(sprintf("  T1 (CS bounds):           %s\n", verdict_T1))
cat(sprintf("  T2 (N-stability):         %d/%d groups\n", n_stable, nrow(t2)))
cat(sprintf("  T3 (M-scaling):           see table; A_hat ~ Var[chi2_com] + 2*sum(lambda^2)\n"))
cat(sprintf("  T4 (imputer-agnosticism): see PMM vs amelia tables\n"))
