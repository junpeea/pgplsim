test_that("Poisson LuGPLSIM fit has expected structure", {

  dat <- make_poisson_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = poisson(),
    offset = dat$exposure,
    M = 4,
    maxit = 30,
    tol = 1e-5
  )

  expect_s3_class(
    fit,
    "LuGPLSIM"
  )

  expect_length(
    fit$alpha_hat,
    ncol(dat$Z)
  )

  expect_length(
    fit$beta_hat,
    ncol(dat$X)
  )

  expect_length(
    fit$fitted_mean,
    length(dat$y)
  )

  expect_true(
    all(is.finite(fit$fitted_mean))
  )

  expect_true(
    all(fit$fitted_mean > 0)
  )
})
