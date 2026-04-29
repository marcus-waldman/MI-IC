#' Meng-Rubin (1992) D_3 Statistic for One Nested Pair
#'
#' Given per-imputation log-likelihoods at own MLEs and at pooled
#' estimates for two nested models (full/reduced), computes the Meng-Rubin
#' D_3 statistic. Returns the scalar r_L correction, the averaged LRT
#' statistics at own MLEs (\code{d_bar}) and at pooled (\code{d_tilde_bar}),
#' the D_3 value itself, and the corrected chi-square
#' \eqn{\chi^2_{D_3} = \bar{\tilde d}/(1 + r_L)}.
#'
#' @param loglik_full_own Per-imputation \eqn{\ell_m(\hat\theta_m^{\text{full}})}.
#' @param loglik_full_pooled Per-imputation \eqn{\ell_m(\bar\theta^{\text{full}})}.
#' @param loglik_red_own Per-imputation \eqn{\ell_m(\hat\theta_m^{\text{reduced}})}.
#' @param loglik_red_pooled Per-imputation \eqn{\ell_m(\bar\theta^{\text{reduced}})}.
#' @param k df difference = npar(full) - npar(reduced).
#' @return List: \code{d_bar}, \code{d_tilde_bar}, \code{r_L},
#'   \code{D3}, \code{chi2_D3}, \code{k}, \code{M}.
compute_D3 <- function(loglik_full_own, loglik_full_pooled,
                       loglik_red_own,  loglik_red_pooled, k) {

  keep <- !is.na(loglik_full_own) & !is.na(loglik_full_pooled) &
          !is.na(loglik_red_own)  & !is.na(loglik_red_pooled)
  loglik_full_own    <- loglik_full_own[keep]
  loglik_full_pooled <- loglik_full_pooled[keep]
  loglik_red_own     <- loglik_red_own[keep]
  loglik_red_pooled  <- loglik_red_pooled[keep]

  M <- length(loglik_full_own)
  if (M < 2 || is.na(k) || k <= 0) {
    return(list(d_bar = NA_real_, d_tilde_bar = NA_real_,
                r_L = NA_real_, D3 = NA_real_, chi2_D3 = NA_real_,
                k = k, M = M))
  }

  d_m       <- 2 * (loglik_full_own    - loglik_red_own)
  d_tilde_m <- 2 * (loglik_full_pooled - loglik_red_pooled)

  d_bar       <- mean(d_m)
  d_tilde_bar <- mean(d_tilde_m)

  r_L <- (M + 1) / (k * (M - 1)) * (d_bar - d_tilde_bar)
  # r_L can be slightly negative due to sampling noise; floor at 0
  r_L <- max(r_L, 0)

  D3      <- d_tilde_bar / (k * (1 + r_L))
  chi2_D3 <- d_tilde_bar / (1 + r_L)

  list(d_bar = d_bar, d_tilde_bar = d_tilde_bar, r_L = r_L,
       D3 = D3, chi2_D3 = chi2_D3, k = k, M = M)
}


#' v4.5 Finite-M Scaled-Shifted Chi-Square Correction
#'
#' Implements the v4.5 §13 closed-form correction for the MI-corrected
#' goodness-of-fit chi-square:
#'   Var[chi2_MI] = Var[chi2_com] + 4 * tr(RIV_perp) + 2 * sum_j lambda_j^2
#'                                                     |_perp + O(1/M)
#' where {lambda_j} are eigenvalues of RIV_perp = I^{imp,mis}_{perp,perp} *
#' (I^{imp,obs}_{perp,perp})^{-1}, restricted to the orthogonal complement of
#' Theta in Phi.
#'
#' The Asparouhov-Muthen-style scaled-shifted correction:
#'   chi2_MI_corr = a * chi2_MI + b
#'   a = sqrt(2 * df / (2 * df + 4 * tr_perp + 2 * sum_lambda_sq_perp))
#'   b = df * (1 - a)
#'
#' By construction, E[chi2_MI_corr] = df and Var[chi2_MI_corr] = 2 df, so
#' chi2_MI_corr can be referenced to a chi^2_df distribution.
#'
#' Approximation used for sum_lambda_sq_perp: leading-order block subtraction
#'   sum_lambda_sq |_perp ~= sum_lambda_sq |_Msat - sum_lambda_sq |_M
#' This neglects the cross-block contribution of order O(1/N); see v4.5 §13
#' for discussion.
#'
#' @param chi2_MI Scalar v4 first-moment-corrected chi-square statistic.
#' @param df Degrees of freedom (= p_Msat - p_M).
#' @param tr_RIV_M Trace of RIV for analysis model.
#' @param tr_RIV_Msat Trace of RIV for saturated reference.
#' @param sum_lambda_sq_M Sum of squared RIV eigenvalues for analysis model.
#' @param sum_lambda_sq_Msat Sum of squared RIV eigenvalues for saturated.
#' @return List with \code{chi2_MI_corr}, \code{a}, \code{b}, \code{tr_perp},
#'   \code{sum_lambda_sq_perp}, \code{var_predicted}.
compute_chi2_MI_corrected <- function(chi2_MI, df,
                                       tr_RIV_M, tr_RIV_Msat,
                                       sum_lambda_sq_M, sum_lambda_sq_Msat) {

  na_result <- list(chi2_MI_corr = NA_real_, a = NA_real_, b = NA_real_,
                    tr_perp = NA_real_, sum_lambda_sq_perp = NA_real_,
                    var_predicted = NA_real_)

  if (is.na(df) || df <= 0 || is.na(chi2_MI)) return(na_result)
  if (any(is.na(c(tr_RIV_M, tr_RIV_Msat,
                   sum_lambda_sq_M, sum_lambda_sq_Msat)))) return(na_result)

  tr_perp           <- tr_RIV_Msat - tr_RIV_M
  sum_lambda_sq_perp <- max(sum_lambda_sq_Msat - sum_lambda_sq_M, 0)

  var_target    <- 2 * df
  var_excess    <- 4 * tr_perp + 2 * sum_lambda_sq_perp
  var_predicted <- var_target + var_excess

  if (var_predicted <= 0) return(na_result)

  a <- sqrt(var_target / var_predicted)
  b <- df * (1 - a)
  chi2_corr <- a * chi2_MI + b

  list(
    chi2_MI_corr       = chi2_corr,
    a                  = a,
    b                  = b,
    tr_perp            = tr_perp,
    sum_lambda_sq_perp = sum_lambda_sq_perp,
    var_predicted      = var_predicted
  )
}


#' Compute Chi-Squares vs Saturated and Finish MR_DEVIANCE
#'
#' For each candidate model \eqn{M_j} (rows of \code{dev_df} other than
#' "Msat"), computes three chi-squares against the saturated reference
#' and fills in the \code{MR_DEVIANCE} column of \code{dev_df}.
#'
#' \itemize{
#'   \item \code{chi2_com} = DEV_com(M_j) - DEV_com(Msat)
#'   \item \code{chi2_MI}  = MI_DEVIANCE(M_j) - MI_DEVIANCE(Msat)
#'   \item \code{chi2_D3}  = Meng-Rubin (1992) corrected chi-square.
#'   \item \code{df}       = npar(Msat) - npar(M_j)
#' }
#'
#' MR_DEVIANCE is anchored at \eqn{-2\bar Q_{\text{sat}}(\bar\theta_{\text{sat}})}
#' for the saturated row; candidate rows get \code{anchor + chi2_D3}.
#'
#' @param complete_fits Named list from \code{\link{fit_complete}}.
#' @param mi_fits Named list from \code{\link{fit_mi_models}}.
#' @param dev_df data.frame from \code{\link{compute_deviances}}.
#' @return List with two elements: \code{chi2_df} (rows = candidate
#'   models, columns chi2_com / chi2_MI / chi2_D3 / df) and
#'   \code{dev_df} (original with MR_DEVIANCE column filled).
compute_chi_squares <- function(complete_fits, mi_fits, dev_df) {

  sat_name <- "Msat"
  if (!(sat_name %in% rownames(dev_df))) {
    stop("compute_chi_squares: dev_df must include a 'Msat' row.")
  }

  # Compute chi-squares for every non-Msat row, including Mnull if present
  # (its chi-square is needed for CFI/TLI denominators). The IC pipeline
  # filters out Mnull separately (compute_all_models_ic operates only on
  # the candidate model list).
  candidates <- setdiff(rownames(dev_df), sat_name)

  sat_mi <- mi_fits[[sat_name]]

  # MR_DEVIANCE anchor at saturated: -2 * mean_m l_m(theta_bar_sat)
  sat_pooled_logliks <- if (isTRUE(sat_mi$success)) {
    sat_mi$logliks_at_pooled[!is.na(sat_mi$logliks_at_pooled)]
  } else {
    numeric(0)
  }
  mr_sat <- if (length(sat_pooled_logliks) > 0) {
    -2 * mean(sat_pooled_logliks)
  } else NA_real_

  dev_df[sat_name, "MR_DEVIANCE"] <- mr_sat

  npar_sat       <- dev_df[sat_name, "npar"]
  dev_com_sat    <- dev_df[sat_name, "DEV_com"]
  dev_adhoc_sat  <- dev_df[sat_name, "DEV_adhoc"]
  mi_dev_sat     <- dev_df[sat_name, "MI_DEVIANCE"]
  tr_RIV_sat     <- if (isTRUE(sat_mi$success)) sat_mi$tr_RIV else NA_real_
  sumlambsq_sat  <- if (isTRUE(sat_mi$success)) sat_mi$sum_lambda_sq else NA_real_

  chi2_rows <- lapply(candidates, function(mname) {
    mf_j   <- mi_fits[[mname]]
    npar_j <- dev_df[mname, "npar"]
    k      <- if (!is.na(npar_sat) && !is.na(npar_j)) npar_sat - npar_j else NA_real_

    chi2_com_j <- if (!is.na(dev_com_sat) && !is.na(dev_df[mname, "DEV_com"])) {
      dev_df[mname, "DEV_com"] - dev_com_sat
    } else NA_real_

    chi2_adhoc_j <- if (!is.na(dev_adhoc_sat) && !is.na(dev_df[mname, "DEV_adhoc"])) {
      dev_df[mname, "DEV_adhoc"] - dev_adhoc_sat
    } else NA_real_

    chi2_MI_j <- if (!is.na(mi_dev_sat) && !is.na(dev_df[mname, "MI_DEVIANCE"])) {
      dev_df[mname, "MI_DEVIANCE"] - mi_dev_sat
    } else NA_real_

    chi2_D3_j <- NA_real_
    if (isTRUE(mf_j$success) && isTRUE(sat_mi$success) &&
        !is.na(k) && k > 0) {
      d3 <- compute_D3(
        loglik_full_own    = sat_mi$logliks,
        loglik_full_pooled = sat_mi$logliks_at_pooled,
        loglik_red_own     = mf_j$logliks,
        loglik_red_pooled  = mf_j$logliks_at_pooled,
        k                  = k
      )
      chi2_D3_j <- d3$chi2_D3
    }

    # v4.5 §13 finite-M scaled-shifted correction
    tr_RIV_j      <- if (isTRUE(mf_j$success)) mf_j$tr_RIV else NA_real_
    sumlambsq_j   <- if (isTRUE(mf_j$success)) mf_j$sum_lambda_sq else NA_real_
    corr_j <- compute_chi2_MI_corrected(
      chi2_MI            = chi2_MI_j,
      df                 = k,
      tr_RIV_M           = tr_RIV_j,
      tr_RIV_Msat        = tr_RIV_sat,
      sum_lambda_sq_M    = sumlambsq_j,
      sum_lambda_sq_Msat = sumlambsq_sat
    )

    data.frame(chi2_com = chi2_com_j, chi2_adhoc = chi2_adhoc_j,
               chi2_MI = chi2_MI_j, chi2_MI_corr = corr_j$chi2_MI_corr,
               chi2_D3 = chi2_D3_j, df = k,
               sum_lambda_sq_M    = sumlambsq_j,
               sum_lambda_sq_Msat = sumlambsq_sat,
               tr_perp            = corr_j$tr_perp,
               sum_lambda_sq_perp = corr_j$sum_lambda_sq_perp,
               var_predicted      = corr_j$var_predicted,
               a_scale            = corr_j$a,
               b_shift            = corr_j$b,
               stringsAsFactors = FALSE)
  })

  chi2_df <- do.call(rbind, chi2_rows)
  rownames(chi2_df) <- candidates

  # Populate MR_DEVIANCE for candidates: anchor + chi2_D3
  dev_df[candidates, "MR_DEVIANCE"] <- mr_sat + chi2_df$chi2_D3

  list(chi2_df = chi2_df, dev_df = dev_df)
}
