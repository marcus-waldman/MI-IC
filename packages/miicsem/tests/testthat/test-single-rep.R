test_that("single replication runs end-to-end with deviances and chi-squares", {
  skip_on_cran()

  config <- get_config(n_reps = 1L, master_seed = 42L)

  models_with_sat <- c(get_sim1_models(),
                       list(Msat = get_saturated_model(config$var_names)))
  pop_starts <- compute_pop_starts(config$sigma_pop, models_with_sat)

  result <- run_one_rep(
    rep_id      = 1,
    n           = 250,
    miss_rate   = 0.25,
    config      = config,
    seed_data   = 42L,
    seed_ampute = 43L,
    seed_impute = 44L,
    pop_starts  = pop_starts
  )

  expect_false(is.null(result))
  expect_true(all(c("ic_df", "dev_df", "chi2_df",
                    "selections", "convergence", "tr_RIVs") %in% names(result)))

  # ic_df: 12 candidates x 7 IC methods
  expect_equal(nrow(result$ic_df), 12L)
  expect_equal(ncol(result$ic_df), 7L)
  expect_true(result$selections["AIC_com"] == "M1")

  # dev_df: 13 rows (12 candidates + Msat) x 7 cols
  expect_equal(nrow(result$dev_df), 13L)
  expect_true("Msat" %in% rownames(result$dev_df))
  expect_true(all(c("DEV_com", "DEV_adhoc", "MI_DEVIANCE", "MR_DEVIANCE",
                    "npar", "tr_RIV", "tr_RIV_fiml") %in%
                    colnames(result$dev_df)))
  # Meanstructure adds 9 intercepts; Msat = 9*10/2 covs + 9 means = 54
  expect_equal(result$dev_df["Msat", "npar"], 54)

  # tr_RIV_fiml populated and positive for at least some models
  expect_true(any(!is.na(result$dev_df$tr_RIV_fiml)))
  non_na <- result$dev_df$tr_RIV_fiml[!is.na(result$dev_df$tr_RIV_fiml)]
  expect_true(all(non_na > -0.5))   # allow slight negative from finite-N

  # chi2_df: 12 candidates x 4 cols
  expect_equal(nrow(result$chi2_df), 12L)
  expect_true(all(c("chi2_com", "chi2_MI", "chi2_D3", "df") %in%
                    colnames(result$chi2_df)))

  # chi-squares non-negative
  expect_true(all(result$chi2_df$chi2_com >= -0.01, na.rm = TRUE))
  expect_true(all(result$chi2_df$chi2_MI  >= -0.01, na.rm = TRUE))
  expect_true(all(result$chi2_df$chi2_D3  >= -0.01, na.rm = TRUE))

  # df = 54 - npar_j for every candidate (saturated npar = 54 with means)
  expect_true(all(
    result$chi2_df$df ==
      54 - result$dev_df[rownames(result$chi2_df), "npar"],
    na.rm = TRUE
  ))

  # MR_DEVIANCE populated for all candidates + saturated
  expect_false(any(is.na(result$dev_df$MR_DEVIANCE)))
})

test_that("get_config returns expected structure", {
  cfg <- get_config()
  expect_equal(dim(cfg$sigma_pop), c(9L, 9L))
  expect_equal(cfg$M, 20L)
  expect_equal(cfg$master_seed, 32897891L)
})

test_that("get_sim1_models returns 12 models", {
  models <- get_sim1_models()
  expect_length(models, 12L)
  expect_equal(names(models), paste0("M", 1:12))
})

test_that("get_saturated_model produces valid lavaan syntax", {
  syn <- get_saturated_model()
  expect_type(syn, "character")
  expect_true(grepl("y1 ~~ y1", syn))
  expect_true(grepl("y9 ~~ y9", syn))
})

test_that("compute_D3 gives correct formula result", {
  # Contrived example: 5 imputations, k = 3
  set.seed(1)
  M <- 5; k <- 3
  l_full_own    <- rnorm(M, -100)
  l_red_own     <- l_full_own - 5       # LRT each imputation ~ 10
  l_full_pooled <- l_full_own - 0.2     # pooled fits slightly worse
  l_red_pooled  <- l_red_own  - 0.2

  res <- compute_D3(l_full_own, l_full_pooled, l_red_own, l_red_pooled, k)

  d_m       <- 2 * (l_full_own    - l_red_own)
  d_tilde_m <- 2 * (l_full_pooled - l_red_pooled)
  expected_d_bar       <- mean(d_m)
  expected_d_tilde_bar <- mean(d_tilde_m)
  expected_r_L <- max((M + 1) / (k * (M - 1)) *
                        (expected_d_bar - expected_d_tilde_bar), 0)
  expected_D3  <- expected_d_tilde_bar / (k * (1 + expected_r_L))

  expect_equal(res$d_bar,       expected_d_bar)
  expect_equal(res$d_tilde_bar, expected_d_tilde_bar)
  expect_equal(res$r_L,         expected_r_L)
  expect_equal(res$D3,          expected_D3)
  expect_equal(res$chi2_D3,     k * expected_D3)
})
