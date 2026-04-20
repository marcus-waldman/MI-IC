# ============================================================================
# Results Analysis for Study 2 (SEM)
# ============================================================================

#' Load and Combine All Condition Results
#'
#' @param results_dir Path to results directory
#' @return data.frame with one row per replication x condition
load_results <- function(results_dir = "simulations/study2-bollen-sem/results") {
  combined_file <- file.path(results_dir, "results_combined.rds")
  if (!file.exists(combined_file)) {
    stop("No combined results file found. Run simulation first.")
  }
  all_results <- readRDS(combined_file)

  rows <- list()
  for (cond_label in names(all_results)) {
    cond_reps <- all_results[[cond_label]]
    for (rep_result in cond_reps) {
      if (is.null(rep_result)) next
      rows[[length(rows) + 1]] <- data.frame(
        n         = rep_result$n,
        miss_rate = rep_result$miss_rate,
        rep_id    = rep_result$rep_id,
        t(rep_result$selections),
        tr_RIV_M1 = rep_result$tr_RIVs["M1"],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


#' Compute Selection Accuracy Table
#'
#' @param results_df data.frame from load_results()
#' @param ic_methods Character vector of IC method column names
#' @return data.frame: condition x method -> % selecting M1
selection_accuracy_table <- function(results_df, ic_methods = NULL) {
  if (is.null(ic_methods)) {
    ic_methods <- c("AIC_com", "BIC_com", "AIC_adhoc", "BIC_adhoc",
                    "AICcd", "MI_AIC", "MI_BIC")
  }

  conditions <- unique(results_df[, c("n", "miss_rate")])
  conditions <- conditions[order(conditions$n, conditions$miss_rate), ]

  acc_rows <- lapply(seq_len(nrow(conditions)), function(i) {
    cond <- conditions[i, ]
    sub <- results_df[results_df$n == cond$n &
                      results_df$miss_rate == cond$miss_rate, ]
    n_valid <- nrow(sub)

    acc <- vapply(ic_methods, function(m) {
      if (!(m %in% colnames(sub))) return(NA_real_)
      sum(sub[[m]] == "M1", na.rm = TRUE) / n_valid * 100
    }, numeric(1))

    data.frame(N = cond$n, miss_rate = cond$miss_rate,
               n_reps = n_valid, t(acc),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, acc_rows)
}


#' Compute Full Selection Frequency Table
#'
#' @param results_df data.frame from load_results()
#' @param ic_method Character. Which IC method to tabulate
#' @return data.frame: condition x model -> selection count
selection_frequency_table <- function(results_df, ic_method = "MI_AIC") {
  model_names <- paste0("M", 1:12)
  conditions <- unique(results_df[, c("n", "miss_rate")])
  conditions <- conditions[order(conditions$n, conditions$miss_rate), ]

  freq_rows <- lapply(seq_len(nrow(conditions)), function(i) {
    cond <- conditions[i, ]
    sub <- results_df[results_df$n == cond$n &
                      results_df$miss_rate == cond$miss_rate, ]
    counts <- table(factor(sub[[ic_method]], levels = model_names))
    data.frame(N = cond$n, miss_rate = cond$miss_rate,
               t(as.numeric(counts)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, freq_rows)
  colnames(out)[3:14] <- model_names
  out
}


#' Summarize tr(RIV) Distributions by Condition
#'
#' @param results_df data.frame from load_results()
#' @return data.frame with mean, median, sd, q05, q95 of tr(RIV) for M1
tr_riv_summary <- function(results_df) {
  conditions <- unique(results_df[, c("n", "miss_rate")])
  conditions <- conditions[order(conditions$n, conditions$miss_rate), ]

  summ_rows <- lapply(seq_len(nrow(conditions)), function(i) {
    cond <- conditions[i, ]
    vals <- results_df$tr_RIV_M1[results_df$n == cond$n &
                                  results_df$miss_rate == cond$miss_rate]
    vals <- vals[!is.na(vals)]
    data.frame(
      N         = cond$n,
      miss_rate = cond$miss_rate,
      mean      = mean(vals),
      median    = median(vals),
      sd        = sd(vals),
      q05       = quantile(vals, 0.05),
      q95       = quantile(vals, 0.95)
    )
  })
  do.call(rbind, summ_rows)
}


#' Plot Selection Accuracy by N
#'
#' @param acc_table data.frame from selection_accuracy_table()
#' @param output_file Path to save PDF (NULL for screen)
plot_accuracy_by_n <- function(acc_table, output_file = NULL) {
  if (!is.null(output_file)) pdf(output_file, width = 10, height = 6)

  ic_methods <- c("AIC_com", "AIC_adhoc", "AICcd", "MI_AIC")
  colors <- c("black", "blue", "red", "darkgreen")
  ltys   <- c(1, 2, 4, 1)

  miss_rates <- sort(unique(acc_table$miss_rate))
  par(mfrow = c(1, length(miss_rates)), mar = c(4, 4, 3, 1))

  for (mr in miss_rates) {
    sub <- acc_table[acc_table$miss_rate == mr, ]

    plot(NULL, xlim = range(sub$N), ylim = c(0, 100),
         xlab = "N", ylab = "% Selecting M1",
         main = sprintf("Miss rate = %.0f%%", mr * 100),
         log = "x")

    for (j in seq_along(ic_methods)) {
      m <- ic_methods[j]
      if (m %in% colnames(sub)) {
        lines(sub$N, sub[[m]], col = colors[j], lty = ltys[j], lwd = 2)
        points(sub$N, sub[[m]], col = colors[j], pch = 16)
      }
    }

    if (mr == miss_rates[1]) {
      legend("bottomright", legend = ic_methods, col = colors,
             lty = ltys, lwd = 2, cex = 0.8)
    }
  }

  if (!is.null(output_file)) dev.off()
}


# --- Main execution ---
if (sys.nframe() == 0) {
  results_df <- load_results()

  cat("=== Selection Accuracy (% selecting M1) ===\n")
  acc <- selection_accuracy_table(results_df)
  print(round(acc, 1))

  cat("\n=== tr(RIV) Summary for M1 ===\n")
  tr_summ <- tr_riv_summary(results_df)
  print(round(tr_summ, 3))

  cat("\n=== Selection Frequency (MI_AIC) ===\n")
  freq <- selection_frequency_table(results_df, "MI_AIC")
  print(freq)

  plot_accuracy_by_n(acc, "simulations/study2-bollen-sem/results/accuracy_plot.pdf")
  cat("\nPlot saved.\n")
}
