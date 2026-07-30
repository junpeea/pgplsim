test_that("binomial LuGPLSIM fit has expected structure", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
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

  expect_true(
    length(fit$gamma_hat) > 0
  )

  expect_true(
    is.numeric(fit$lambda)
  )

  expect_length(
    fit$fitted_mean,
    length(dat$y)
  )

  expect_length(
    fit$u_hat,
    length(dat$y)
  )
})

test_that("binomial fitted probabilities are valid", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  expect_true(
    all(is.finite(fit$fitted_mean))
  )

  expect_true(
    all(fit$fitted_mean >= 0)
  )

  expect_true(
    all(fit$fitted_mean <= 1)
  )
})
