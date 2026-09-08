test_that("model-based covariance has correct dimensions", {

  dat <- make_binomial_test_data()

  fit <- pgplsim:::pgplsim_fit(
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

  fit <- pgplsim:::pgplsim_fit(
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

  fit <- pgplsim:::pgplsim_fit(
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


test_that("finite-sample covariance is available", {

  dat <- make_binomial_test_data(n = 120, seed = 1401)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  V <- vcov(
    fit,
    type = "finite",
    component = "parameters"
  )

  expected_dimension <- ncol(dat$Z) + ncol(dat$X)

  expect_equal(dim(V), c(expected_dimension, expected_dimension))
  expect_equal(V, t(V), tolerance = 1e-8)
  expect_true(all(is.finite(V)))
})


test_that("leverage-adjusted covariance is available", {

  dat <- make_binomial_test_data(n = 120, seed = 1402)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  V <- vcov(
    fit,
    type = "leverage",
    component = "parameters"
  )

  expected_dimension <- ncol(dat$Z) + ncol(dat$X)

  expect_equal(dim(V), c(expected_dimension, expected_dimension))
  expect_equal(V, t(V), tolerance = 1e-8)
  expect_true(all(is.finite(V)))
})


test_that("bootstrap covariance can be extracted through vcov", {

  dat <- make_binomial_test_data(n = 100, seed = 1403)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot <- bootstrap_pgplsim(
    object = fit,
    B = 3,
    method = "stratified",
    seed = 2026,
    max_attempts = 12,
    trace = FALSE
  )

  V <- vcov(
    fit,
    type = "bootstrap",
    component = "parameters",
    bootstrap = boot
  )

  expected_dimension <- ncol(dat$Z) + ncol(dat$X)

  expect_equal(dim(V), c(expected_dimension, expected_dimension))
  expect_equal(V, t(V), tolerance = 1e-8)
  expect_true(all(is.finite(V)))
})


test_that("bootstrap covariance requires a bootstrap object", {

  dat <- make_binomial_test_data(n = 80, seed = 1404)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  expect_error(
    vcov(fit, type = "bootstrap"),
    "bootstrap result is required"
  )
})


test_that("unsupported covariance types are rejected", {

  dat <- make_binomial_test_data(n = 80, seed = 1405)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  expect_error(
    vcov(fit, type = "not-a-method"),
    "should be one of"
  )
})
