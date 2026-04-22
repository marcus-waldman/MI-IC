# ============================================================================
# Rebuild IC selections under alternative correction formulas
# ============================================================================
# Uses stored per-rep dev_df (MI_DEVIANCE, tr_RIV, npar, DEV_adhoc, DEV_com)
# to compute NEW selection results without re-running the simulation.
#
# Variants computed (on the deviance scale; all add the usual 2p or p*log(n)):
#
#   AIC_com       = DEV_com             + 2p                 (oracle)
#   AIC_adhoc     = DEV_adhoc           + 2p                 (mean at each own MLE)
#   AICcd         = DEV_adhoc           + 2p + 2*tr_RIV      (Cavanaugh-Shumway)
#   MI_AIC_add    = -2Q_bar             + 2p +   tr_RIV      (current miicsem)
#   MI_AIC_sub    = -2Q_bar             + 2p -   tr_RIV      (sign-flipped hypothesis)
#   MI_AIC_nocor  = -2Q_bar             + 2p                 (no RIV correction)
#   MI_AIC_sub2   = -2Q_bar             + 2p - 2*tr_RIV      (double-subtract)
#
# BIC analogs replace 2p with p*log(n) and use the same corrections.
#
# Selection rule: arg min over 12 candidate models. Accuracy = % picking M1.
# ============================================================================

RESULTS <- "hpc/results-full-M100"
FIGURES <- "hpc/figures-M100"

cat("Loading combined results...\n")
all_res <- readRDS(file.path(RESULTS, "results_combined.rds"))
cat(sprintf("  %d condition lists\n", length(all_res)))

IC_NAMES <- c(
  "AIC_com", "BIC_com",
  "AIC_adhoc", "BIC_adhoc",
  "AICcd",
  "MI_AIC_add",    "MI_BIC_add",
  "MI_AIC_sub",    "MI_BIC_sub",
  "MI_AIC_nocor",  "MI_BIC_nocor",
  "MI_AIC_sub2",   "MI_BIC_sub2"
)

# ---- helper: compute all ICs for a single rep's dev_df --------------------
compute_variants <- function(dev_df_rep, n) {
  # keep only candidate models M1..M12, in stable order
  cand <- paste0("M", 1:12)
  d <- dev_df_rep[cand, , drop = FALSE]

  raw_pool <- d$MI_DEVIANCE - d$tr_RIV   # -2*Q_bar
  p        <- d$npar
  logn     <- log(n)
  tR       <- d$tr_RIV

  out <- data.frame(
    model      = cand,
    AIC_com       = d$DEV_com   + 2 * p,
    BIC_com       = d$DEV_com   + p * logn,
    AIC_adhoc     = d$DEV_adhoc + 2 * p,
    BIC_adhoc     = d$DEV_adhoc + p * logn,
    AICcd         = d$DEV_adhoc + 2 * p + 2 * tR,
    MI_AIC_add    = raw_pool    + 2 * p +     tR,
    MI_BIC_add    = raw_pool    + p * logn + tR,
    MI_AIC_sub    = raw_pool    + 2 * p -     tR,
    MI_BIC_sub    = raw_pool    + p * logn - tR,
    MI_AIC_nocor  = raw_pool    + 2 * p,
    MI_BIC_nocor  = raw_pool    + p * logn,
    MI_AIC_sub2   = raw_pool    + 2 * p - 2 * tR,
    MI_BIC_sub2   = raw_pool    + p * logn - 2 * tR,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

# ---- iterate: build a selection table -------------------------------------
cat("\nRebuilding selections under each variant...\n")
rows <- list()
for (cond_label in names(all_res)) {
  reps <- all_res[[cond_label]]
  for (rr in reps) {
    if (is.null(rr) || isTRUE(rr$failed) || is.null(rr$dev_df)) next
    d <- rr$dev_df
    # Skip if required columns missing / bad
    if (!all(c("DEV_com", "DEV_adhoc", "MI_DEVIANCE",
               "npar", "tr_RIV") %in% colnames(d))) next
    ic <- compute_variants(d, rr$n)

    sel <- vapply(IC_NAMES, function(m) {
      v <- ic[[m]]
      if (all(is.na(v))) return(NA_character_)
      ic$model[which.min(v)]
    }, character(1))

    rows[[length(rows) + 1]] <- data.frame(
      n         = rr$n,
      miss_rate = rr$miss_rate,
      rep_id    = rr$rep_id,
      t(sel),
      stringsAsFactors = FALSE
    )
  }
}
sel_df <- do.call(rbind, rows)
cat(sprintf("  %d reps processed\n", nrow(sel_df)))

# ---- selection accuracy tables -------------------------------------------
acc_table <- function(df, ic_names) {
  conds <- unique(df[, c("n", "miss_rate")])
  conds <- conds[order(conds$n, conds$miss_rate), ]
  out <- lapply(seq_len(nrow(conds)), function(i) {
    cc  <- conds[i, ]
    sub <- df[df$n == cc$n & df$miss_rate == cc$miss_rate, ]
    acc <- vapply(ic_names, function(m) {
      if (!(m %in% colnames(sub))) return(NA_real_)
      sum(sub[[m]] == "M1", na.rm = TRUE) / nrow(sub) * 100
    }, numeric(1))
    data.frame(N = cc$n, miss_rate = cc$miss_rate,
               n_reps = nrow(sub), t(acc),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

acc <- acc_table(sel_df, IC_NAMES)

# ---- AIC-class comparison -------------------------------------------------
cat("\n=========================================================\n")
cat("  AIC-class: % selecting true model M1 (by condition)\n")
cat("=========================================================\n")
aic_cols <- c("N", "miss_rate", "n_reps",
              "AIC_com", "AIC_adhoc", "AICcd",
              "MI_AIC_add", "MI_AIC_sub",
              "MI_AIC_nocor", "MI_AIC_sub2")
print(round(acc[, aic_cols], 1), row.names = FALSE)

# ---- BIC-class comparison -------------------------------------------------
cat("\n=========================================================\n")
cat("  BIC-class: % selecting true model M1 (by condition)\n")
cat("=========================================================\n")
bic_cols <- c("N", "miss_rate", "n_reps",
              "BIC_com", "BIC_adhoc",
              "MI_BIC_add", "MI_BIC_sub",
              "MI_BIC_nocor", "MI_BIC_sub2")
print(round(acc[, bic_cols], 1), row.names = FALSE)

# ---- Distance-to-oracle summary ------------------------------------------
cat("\n=========================================================\n")
cat("  Distance from oracle (|method - AIC_com|) averaged across 12 conditions\n")
cat("=========================================================\n")
mean_abs_gap <- function(x, ref) mean(abs(x - ref))
aic_gap <- data.frame(
  method = c("AIC_adhoc", "AICcd",
             "MI_AIC_add", "MI_AIC_sub",
             "MI_AIC_nocor", "MI_AIC_sub2"),
  mean_abs_gap = c(
    mean_abs_gap(acc$AIC_adhoc,    acc$AIC_com),
    mean_abs_gap(acc$AICcd,        acc$AIC_com),
    mean_abs_gap(acc$MI_AIC_add,   acc$AIC_com),
    mean_abs_gap(acc$MI_AIC_sub,   acc$AIC_com),
    mean_abs_gap(acc$MI_AIC_nocor, acc$AIC_com),
    mean_abs_gap(acc$MI_AIC_sub2,  acc$AIC_com)
  )
)
aic_gap <- aic_gap[order(aic_gap$mean_abs_gap), ]
print(aic_gap, row.names = FALSE, digits = 3)

bic_gap <- data.frame(
  method = c("BIC_adhoc",
             "MI_BIC_add", "MI_BIC_sub",
             "MI_BIC_nocor", "MI_BIC_sub2"),
  mean_abs_gap = c(
    mean_abs_gap(acc$BIC_adhoc,    acc$BIC_com),
    mean_abs_gap(acc$MI_BIC_add,   acc$BIC_com),
    mean_abs_gap(acc$MI_BIC_sub,   acc$BIC_com),
    mean_abs_gap(acc$MI_BIC_nocor, acc$BIC_com),
    mean_abs_gap(acc$MI_BIC_sub2,  acc$BIC_com)
  )
)
bic_gap <- bic_gap[order(bic_gap$mean_abs_gap), ]
cat("\nBIC-class gap to BIC_com:\n")
print(bic_gap, row.names = FALSE, digits = 3)

# ---- Deviance-level bias under each correction (M1 only) ------------------
# Sanity: check that the subtracted-correction renders MI_DEVIANCE unbiased
# for DEV_com at the true model.
cat("\n=========================================================\n")
cat("  Deviance bias (M1 only) under each pooled-deviance definition\n")
cat("    bias = E[def(-2Q_bar + k*tr_RIV)] - E[DEV_com]\n")
cat("=========================================================\n")

bias_rows <- list()
for (cond_label in names(all_res)) {
  for (rr in all_res[[cond_label]]) {
    if (is.null(rr) || is.null(rr$dev_df)) next
    d <- rr$dev_df
    if (!"M1" %in% rownames(d)) next
    r <- d["M1", ]
    if (any(is.na(c(r$DEV_com, r$MI_DEVIANCE, r$tr_RIV)))) next
    rawpool <- r$MI_DEVIANCE - r$tr_RIV
    bias_rows[[length(bias_rows) + 1]] <- data.frame(
      n         = rr$n,
      miss_rate = rr$miss_rate,
      tr_RIV    = r$tr_RIV,
      d_raw     = rawpool                - r$DEV_com,
      d_add     = (rawpool +     r$tr_RIV) - r$DEV_com,  # MI_DEV (current)
      d_sub     = (rawpool -     r$tr_RIV) - r$DEV_com,  # hypothesis
      d_sub2    = (rawpool - 2 * r$tr_RIV) - r$DEV_com,
      stringsAsFactors = FALSE
    )
  }
}
bdf <- do.call(rbind, bias_rows)
bias_summ <- aggregate(
  cbind(tr_RIV, d_raw, d_add, d_sub, d_sub2) ~ n + miss_rate,
  data = bdf, FUN = mean
)
bias_summ <- bias_summ[order(bias_summ$n, bias_summ$miss_rate), ]
print(round(bias_summ, 3), row.names = FALSE)

# ---- save -----------------------------------------------------------------
write.csv(acc,        file.path(FIGURES, "accuracy_hypothesized_corrections.csv"), row.names = FALSE)
write.csv(aic_gap,    file.path(FIGURES, "aic_gap_to_oracle.csv"),                  row.names = FALSE)
write.csv(bic_gap,    file.path(FIGURES, "bic_gap_to_oracle.csv"),                  row.names = FALSE)
write.csv(bias_summ,  file.path(FIGURES, "dev_bias_by_correction_M1.csv"),          row.names = FALSE)

cat(sprintf("\n[saved] %s/accuracy_hypothesized_corrections.csv\n", FIGURES))
cat(sprintf("[saved] %s/aic_gap_to_oracle.csv\n", FIGURES))
cat(sprintf("[saved] %s/bic_gap_to_oracle.csv\n", FIGURES))
cat(sprintf("[saved] %s/dev_bias_by_correction_M1.csv\n", FIGURES))
cat("\nDone.\n")
