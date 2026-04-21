#' Load and Combine Selection Results
#'
#' Reads \code{results_combined.rds} produced by \code{\link{run_simulation}}
#' and returns a flat data.frame with one row per rep x condition —
#' focused on model selection (which model won under each IC).
#'
#' @param results_dir Path to results directory.
#' @return A data.frame with columns n, miss_rate, rep_id, one column per
#'   IC method (selected model name), and tr_RIV_M1.
#' @export
load_results <- function(results_dir) {
  combined_file <- file.path(results_dir, "results_combined.rds")
  if (!file.exists(combined_file)) {
    stop("No combined results file found at: ", combined_file)
  }
  all_results <- readRDS(combined_file)

  rows <- list()
  for (cond_label in names(all_results)) {
    cond_reps <- all_results[[cond_label]]
    for (rep_result in cond_reps) {
      if (is.null(rep_result) || isTRUE(rep_result$failed)) next
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


#' Load Long-Form Deviance Table
#'
#' Returns one row per (rep x condition x model) with the four deviance
#' quantities plus npar and tr_RIV.
#'
#' @param results_dir Path to results directory.
#' @return A data.frame with columns n, miss_rate, rep_id, model,
#'   DEV_com, DEV_adhoc, MI_DEVIANCE, MR_DEVIANCE, npar, tr_RIV.
#' @export
load_deviances <- function(results_dir) {
  combined_file <- file.path(results_dir, "results_combined.rds")
  if (!file.exists(combined_file)) {
    stop("No combined results file found at: ", combined_file)
  }
  all_results <- readRDS(combined_file)

  rows <- list()
  for (cond_label in names(all_results)) {
    for (rep_result in all_results[[cond_label]]) {
      if (is.null(rep_result) || is.null(rep_result$dev_df)) next
      dev_df <- rep_result$dev_df
      rows[[length(rows) + 1]] <- data.frame(
        n         = rep_result$n,
        miss_rate = rep_result$miss_rate,
        rep_id    = rep_result$rep_id,
        model     = rownames(dev_df),
        dev_df,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}


#' Load Long-Form Chi-Square Table
#'
#' Returns one row per (rep x condition x candidate model) with the three
#' chi-square variants plus df.
#'
#' @param results_dir Path to results directory.
#' @return A data.frame with columns n, miss_rate, rep_id, model,
#'   chi2_com, chi2_MI, chi2_D3, df.
#' @export
load_chi_squares <- function(results_dir) {
  combined_file <- file.path(results_dir, "results_combined.rds")
  if (!file.exists(combined_file)) {
    stop("No combined results file found at: ", combined_file)
  }
  all_results <- readRDS(combined_file)

  rows <- list()
  for (cond_label in names(all_results)) {
    for (rep_result in all_results[[cond_label]]) {
      if (is.null(rep_result) || is.null(rep_result$chi2_df)) next
      chi2_df <- rep_result$chi2_df
      rows[[length(rows) + 1]] <- data.frame(
        n         = rep_result$n,
        miss_rate = rep_result$miss_rate,
        rep_id    = rep_result$rep_id,
        model     = rownames(chi2_df),
        chi2_df,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}


#' Summarize Chi-Squares by Condition and Model
#'
#' Mean of each chi-square variant and the oracle \eqn{\chi^2_{\text{com}}}
#' for each (N, miss_rate, model) cell.
#'
#' @param chi2_df Long-form data.frame from \code{\link{load_chi_squares}}.
#' @return data.frame with mean chi2_com, chi2_MI, chi2_D3 per cell.
#' @export
chi2_summary_table <- function(chi2_df) {
  agg <- stats::aggregate(
    cbind(chi2_com, chi2_MI, chi2_D3) ~ n + miss_rate + model,
    data = chi2_df,
    FUN  = function(x) mean(x, na.rm = TRUE)
  )
  agg[order(agg$n, agg$miss_rate, agg$model), ]
}


#' Compute Selection Accuracy Table
#'
#' Percentage of replications in which each IC method selected the
#' true model (M1).
#'
#' @param results_df data.frame from \code{\link{load_results}}.
#' @param ic_methods Character vector of IC method column names.
#' @return data.frame: one row per condition, columns for each method.
#' @export
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


#' Full Selection Frequency Table for a Single IC Method
#'
#' @param results_df data.frame from \code{\link{load_results}}.
#' @param ic_method Character. Which IC method to tabulate.
#' @return data.frame with counts for each of the 12 models.
#' @export
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


#' Summarize tr(RIV) Distributions for M1 by Condition
#'
#' @param results_df data.frame from \code{\link{load_results}}.
#' @return data.frame with mean, median, sd, q05, q95 of tr(RIV).
#' @export
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
      median    = stats::median(vals),
      sd        = stats::sd(vals),
      q05       = stats::quantile(vals, 0.05),
      q95       = stats::quantile(vals, 0.95)
    )
  })
  do.call(rbind, summ_rows)
}


#' Plot Selection Accuracy by Sample Size
#'
#' @param acc_table data.frame from \code{\link{selection_accuracy_table}}.
#' @param output_file Path to save PDF (NULL = screen).
#' @export
plot_accuracy_by_n <- function(acc_table, output_file = NULL) {
  if (!is.null(output_file)) grDevices::pdf(output_file, width = 10, height = 6)

  ic_methods <- c("AIC_com", "AIC_adhoc", "AICcd", "MI_AIC")
  colors     <- c("black", "blue", "red", "darkgreen")
  ltys       <- c(1, 2, 4, 1)

  miss_rates <- sort(unique(acc_table$miss_rate))
  graphics::par(mfrow = c(1, length(miss_rates)), mar = c(4, 4, 3, 1))

  for (mr in miss_rates) {
    sub <- acc_table[acc_table$miss_rate == mr, ]

    graphics::plot(NULL, xlim = range(sub$N), ylim = c(0, 100),
         xlab = "N", ylab = "% Selecting M1",
         main = sprintf("Miss rate = %.0f%%", mr * 100),
         log = "x")

    for (j in seq_along(ic_methods)) {
      m <- ic_methods[j]
      if (m %in% colnames(sub)) {
        graphics::lines(sub$N, sub[[m]], col = colors[j], lty = ltys[j], lwd = 2)
        graphics::points(sub$N, sub[[m]], col = colors[j], pch = 16)
      }
    }

    if (mr == miss_rates[1]) {
      graphics::legend("bottomright", legend = ic_methods, col = colors,
             lty = ltys, lwd = 2, cex = 0.8)
    }
  }

  if (!is.null(output_file)) grDevices::dev.off()
}
