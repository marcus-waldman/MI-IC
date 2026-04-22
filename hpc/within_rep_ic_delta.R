# ============================================================================
# Within-rep IC differences: how faithfully does each correction replicate
# the complete-data IC *gaps* between candidate models?
# ============================================================================
# Per rep, compute Delta IC(M_j, M1) = IC(M_j) - IC(M1) for j = 2..12.
# For each IC method, compare Delta to the oracle Delta_AIC_com (or Delta_BIC_com).
#
# Key diagnostics (per (N, miss_rate, comparison_model)):
#   bias        = E[ Delta_method - Delta_oracle ]
#   rmse        = sqrt(E[ (Delta_method - Delta_oracle)^2 ])
#   sign_agree  = P( sign(Delta_method) == sign(Delta_oracle) )
#   rank_agree  = within-rep winner of method == within-rep winner of oracle
#
# Variants (deviance + penalty, selection = argmin over 12 candidates):
#   MI_AIC_add   = -2Q_bar + 2p +   tr_RIV
#   MI_AIC_sub   = -2Q_bar + 2p -   tr_RIV
#   MI_AIC_nocor = -2Q_bar + 2p
#   AIC_adhoc    = DEV_adhoc + 2p
#   AICcd        = DEV_adhoc + 2p + 2*tr_RIV
# (+ BIC analogs)
# ============================================================================

RESULTS <- "hpc/results-full-M100"
FIGURES <- "hpc/figures-M100"

cat("Loading combined results...\n")
all_res <- readRDS(file.path(RESULTS, "results_combined.rds"))

CANDS <- paste0("M", 1:12)

# ---- per-rep IC matrix for all variants -----------------------------------
ic_matrix <- function(d, n) {
  d <- d[CANDS, , drop = FALSE]
  raw_pool <- d$MI_DEVIANCE - d$tr_RIV   # -2 Q_bar
  p    <- d$npar
  logn <- log(n)
  tR   <- d$tr_RIV

  data.frame(
    AIC_com      = d$DEV_com   + 2 * p,
    BIC_com      = d$DEV_com   + p * logn,
    AIC_adhoc    = d$DEV_adhoc + 2 * p,
    BIC_adhoc    = d$DEV_adhoc + p * logn,
    AICcd        = d$DEV_adhoc + 2 * p + 2 * tR,
    MI_AIC_add   = raw_pool    + 2 * p +     tR,
    MI_AIC_sub   = raw_pool    + 2 * p -     tR,
    MI_AIC_nocor = raw_pool    + 2 * p,
    MI_BIC_add   = raw_pool    + p * logn + tR,
    MI_BIC_sub   = raw_pool    + p * logn - tR,
    MI_BIC_nocor = raw_pool    + p * logn,
    row.names = CANDS,
    stringsAsFactors = FALSE
  )
}

# ---- build long-format Delta IC data --------------------------------------
cat("Building within-rep Delta IC table (this takes a minute)...\n")
rows <- vector("list", length(all_res) * 2000)
k <- 0
for (cond_label in names(all_res)) {
  for (rr in all_res[[cond_label]]) {
    if (is.null(rr) || is.null(rr$dev_df)) next
    d <- rr$dev_df
    if (!all(CANDS %in% rownames(d))) next
    if (any(is.na(d[CANDS, c("DEV_com", "DEV_adhoc",
                             "MI_DEVIANCE", "npar", "tr_RIV")]))) next
    ic <- ic_matrix(d, rr$n)
    # Delta_j_1 for each IC column:  ic[j,] - ic["M1",]
    delta <- sweep(as.matrix(ic), 2, as.numeric(ic["M1", ]), "-")
    # drop M1 row (Delta = 0 by construction)
    delta <- delta[CANDS[-1], , drop = FALSE]
    for (j in seq_len(nrow(delta))) {
      k <- k + 1
      rows[[k]] <- data.frame(
        n         = rr$n,
        miss_rate = rr$miss_rate,
        rep_id    = rr$rep_id,
        alt_model = rownames(delta)[j],
        t(delta[j, ]),
        stringsAsFactors = FALSE
      )
    }
  }
}
rows <- rows[seq_len(k)]
dlong <- do.call(rbind, rows)
rm(rows); gc()
cat(sprintf("  %d rep x alt-model rows\n", nrow(dlong)))

# ---- AIC-side analysis ----------------------------------------------------
aic_methods <- c("AIC_adhoc", "AICcd", "MI_AIC_add", "MI_AIC_sub", "MI_AIC_nocor")

cat("\n=========================================================\n")
cat("  AIC-class: within-rep replication of oracle Delta AIC_com\n")
cat("  Averaged over 11 alt models x reps per condition\n")
cat("=========================================================\n")

aic_rows <- list()
conds <- unique(dlong[, c("n", "miss_rate")])
conds <- conds[order(conds$n, conds$miss_rate), ]
for (i in seq_len(nrow(conds))) {
  cc  <- conds[i, ]
  sub <- dlong[dlong$n == cc$n & dlong$miss_rate == cc$miss_rate, ]
  oracle <- sub$AIC_com
  r <- data.frame(N = cc$n, miss_rate = cc$miss_rate,
                  n_pairs = nrow(sub), stringsAsFactors = FALSE)
  for (m in aic_methods) {
    diff <- sub[[m]] - oracle
    r[[paste0(m, "__bias")]] <- mean(diff)
    r[[paste0(m, "__rmse")]] <- sqrt(mean(diff^2))
    r[[paste0(m, "__sgn")]]  <- mean(sign(sub[[m]]) == sign(oracle)) * 100
  }
  aic_rows[[i]] <- r
}
aic_summ <- do.call(rbind, aic_rows)

# Pretty-print in three blocks: bias, rmse, sign-agreement
print_block <- function(df, suffix, label, digits = 2) {
  cols <- grep(paste0("__", suffix, "$"), colnames(df), value = TRUE)
  tbl  <- df[, c("N", "miss_rate", "n_pairs", cols)]
  colnames(tbl) <- c("N", "miss_rate", "n_pairs",
                     sub(paste0("__", suffix), "", cols))
  cat(sprintf("\n-- %s --\n", label))
  print(round(tbl[, -3], digits), row.names = FALSE)
}

print_block(aic_summ, "bias", "Bias: E[Delta_method - Delta_AIC_com]", 3)
print_block(aic_summ, "rmse", "RMSE: sqrt(E[(Delta_method - Delta_AIC_com)^2])", 2)
print_block(aic_summ, "sgn",  "Sign agreement %: sign(Delta_method) == sign(Delta_AIC_com)", 1)

# ---- BIC-side analysis ----------------------------------------------------
bic_methods <- c("BIC_adhoc", "MI_BIC_add", "MI_BIC_sub", "MI_BIC_nocor")

cat("\n=========================================================\n")
cat("  BIC-class: within-rep replication of oracle Delta BIC_com\n")
cat("=========================================================\n")

bic_rows <- list()
for (i in seq_len(nrow(conds))) {
  cc  <- conds[i, ]
  sub <- dlong[dlong$n == cc$n & dlong$miss_rate == cc$miss_rate, ]
  oracle <- sub$BIC_com
  r <- data.frame(N = cc$n, miss_rate = cc$miss_rate,
                  n_pairs = nrow(sub), stringsAsFactors = FALSE)
  for (m in bic_methods) {
    diff <- sub[[m]] - oracle
    r[[paste0(m, "__bias")]] <- mean(diff)
    r[[paste0(m, "__rmse")]] <- sqrt(mean(diff^2))
    r[[paste0(m, "__sgn")]]  <- mean(sign(sub[[m]]) == sign(oracle)) * 100
  }
  bic_rows[[i]] <- r
}
bic_summ <- do.call(rbind, bic_rows)

print_block(bic_summ, "bias", "Bias: E[Delta_method - Delta_BIC_com]",  3)
print_block(bic_summ, "rmse", "RMSE",                                    2)
print_block(bic_summ, "sgn",  "Sign agreement %",                        1)

# ---- Rank-level: winner agreement within rep ------------------------------
cat("\n=========================================================\n")
cat("  Winner-agreement within rep: P( argmin_method == argmin_oracle )\n")
cat("  (Does the method pick the SAME model the oracle picks in this rep?)\n")
cat("=========================================================\n")

winner_rows <- list()
k <- 0
for (cond_label in names(all_res)) {
  for (rr in all_res[[cond_label]]) {
    if (is.null(rr) || is.null(rr$dev_df)) next
    d <- rr$dev_df
    if (!all(CANDS %in% rownames(d))) next
    if (any(is.na(d[CANDS, c("DEV_com", "DEV_adhoc",
                             "MI_DEVIANCE", "npar", "tr_RIV")]))) next
    ic <- ic_matrix(d, rr$n)

    aic_com_win <- CANDS[which.min(ic$AIC_com)]
    bic_com_win <- CANDS[which.min(ic$BIC_com)]

    k <- k + 1
    winner_rows[[k]] <- data.frame(
      n          = rr$n,
      miss_rate  = rr$miss_rate,
      rep_id     = rr$rep_id,
      win_AIC_com        = aic_com_win,
      win_BIC_com        = bic_com_win,
      win_AIC_adhoc      = CANDS[which.min(ic$AIC_adhoc)],
      win_AICcd          = CANDS[which.min(ic$AICcd)],
      win_MI_AIC_add     = CANDS[which.min(ic$MI_AIC_add)],
      win_MI_AIC_sub     = CANDS[which.min(ic$MI_AIC_sub)],
      win_MI_AIC_nocor   = CANDS[which.min(ic$MI_AIC_nocor)],
      win_BIC_adhoc      = CANDS[which.min(ic$BIC_adhoc)],
      win_MI_BIC_add     = CANDS[which.min(ic$MI_BIC_add)],
      win_MI_BIC_sub     = CANDS[which.min(ic$MI_BIC_sub)],
      win_MI_BIC_nocor   = CANDS[which.min(ic$MI_BIC_nocor)],
      stringsAsFactors = FALSE
    )
  }
}
win_df <- do.call(rbind, winner_rows)

agree_table <- function(df, ref_col, method_cols) {
  conds <- unique(df[, c("n", "miss_rate")])
  conds <- conds[order(conds$n, conds$miss_rate), ]
  do.call(rbind, lapply(seq_len(nrow(conds)), function(i) {
    cc <- conds[i, ]
    sub <- df[df$n == cc$n & df$miss_rate == cc$miss_rate, ]
    pct <- vapply(method_cols, function(m) {
      mean(sub[[m]] == sub[[ref_col]]) * 100
    }, numeric(1))
    data.frame(N = cc$n, miss_rate = cc$miss_rate,
               n_reps = nrow(sub), t(pct),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

aic_agree <- agree_table(win_df, "win_AIC_com",
                         c("win_AIC_adhoc", "win_AICcd",
                           "win_MI_AIC_add", "win_MI_AIC_sub",
                           "win_MI_AIC_nocor"))
cat("\n-- AIC-side: % within-rep winner matches AIC_com --\n")
print(round(aic_agree, 1), row.names = FALSE)

bic_agree <- agree_table(win_df, "win_BIC_com",
                         c("win_BIC_adhoc",
                           "win_MI_BIC_add", "win_MI_BIC_sub",
                           "win_MI_BIC_nocor"))
cat("\n-- BIC-side: % within-rep winner matches BIC_com --\n")
print(round(bic_agree, 1), row.names = FALSE)

# ---- save ----------------------------------------------------------------
write.csv(aic_summ,  file.path(FIGURES, "delta_ic_aic_summary.csv"),  row.names = FALSE)
write.csv(bic_summ,  file.path(FIGURES, "delta_ic_bic_summary.csv"),  row.names = FALSE)
write.csv(aic_agree, file.path(FIGURES, "within_rep_aic_winner_agreement.csv"), row.names = FALSE)
write.csv(bic_agree, file.path(FIGURES, "within_rep_bic_winner_agreement.csv"), row.names = FALSE)

cat(sprintf("\n[saved] %s/delta_ic_{aic,bic}_summary.csv\n", FIGURES))
cat(sprintf("[saved] %s/within_rep_{aic,bic}_winner_agreement.csv\n", FIGURES))
cat("\nDone.\n")
