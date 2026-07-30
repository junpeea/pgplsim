test_that("bootstrap returns expected estimates", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 789
  )

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot <- bootstrap_LuGPLSIM(
    object = fit,
    B = 3,
    method = "stratified",
    seed = 2026,
    max_attempts = 10,
    trace = FALSE
  )

  expect_s3_class(
    boot,
    "bootstrap.LuGPLSIM"
  )

  expect_equal(
    nrow(boot$estimates_alpha),
    3
  )

  expect_equal(
    ncol(boot$estimates_alpha),
    ncol(dat$Z)
  )

  expect_equal(
    nrow(boot$estimates_beta),
    3
  )

  expect_equal(
    ncol(boot$estimates_beta),
    ncol(dat$X)
  )
})
