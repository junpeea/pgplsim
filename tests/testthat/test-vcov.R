test_that("model-based covariance has correct dimensions", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  V_parameters <- vcov(
    fit,
    type = "model",
    component = "parameters"
  )

  expected_dimension <-
    ncol(dat$Z) +
    ncol(dat$X)

  expect_equal(
    nrow(V_parameters),
    expected_dimension
  )

  expect_equal(
    ncol(V_parameters),
    expected_dimension
  )

  expect_equal(
    V_parameters,
    t(V_parameters),
    tolerance = 1e-8
  )
})

test_that("alpha and beta covariance blocks have correct dimensions", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  V_alpha <- vcov(
    fit,
    type = "model",
    component = "alpha"
  )

  V_beta <- vcov(
    fit,
    type = "model",
    component = "beta"
  )

  expect_equal(
    dim(V_alpha),
    c(ncol(dat$Z), ncol(dat$Z))
  )

  expect_equal(
    dim(V_beta),
    c(ncol(dat$X), ncol(dat$X))
  )
})

test_that("sandwich covariance is available", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  V <- vcov(
    fit,
    type = "sandwich",
    component = "parameters"
  )

  expected_dimension <-
    ncol(dat$Z) +
    ncol(dat$X)

  expect_equal(
    dim(V),
    c(expected_dimension, expected_dimension)
  )

  expect_true(
    all(is.finite(V))
  )
})
