test_that("single replication runs end-to-end", {
  skip_on_cran()

  config <- get_config(n_reps = 1L, master_seed = 42L)

  result <- run_one_rep(
    rep_id      = 1,
    n           = 250,
    miss_rate   = 0.25,
    config      = config,
    seed_data   = 42L,
    seed_ampute = 43L,
    seed_impute = 44L
  )

  expect_false(is.null(result))
  expect_true(all(c("ic_df", "selections", "convergence", "tr_RIVs") %in% names(result)))
  expect_equal(nrow(result$ic_df), 12L)
  expect_equal(ncol(result$ic_df), 7L)
  expect_true(result$selections["AIC_com"] == "M1")
  expect_true(all(result$tr_RIVs > 0, na.rm = TRUE))
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
