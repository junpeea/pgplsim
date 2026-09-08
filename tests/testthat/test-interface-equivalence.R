test_that("formula and matrix interfaces give equivalent binomial fits", {

  set.seed(2026)

  n <- 200

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  eta <- 0.4 +
    0.6 * dat$x1 -
    0.3 * dat$x2 +
    0.5 * sin(dat$z1 + dat$z2)

  dat$y <- rbinom(
    n,
    size = 1,
    prob = plogis(eta)
  )

  ctrl <- pgplsim_control(
    maxit = 30,
    tol = 1e-5,
    trace = FALSE
  )

  fit_public <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4,
    lambda = 1,
    control = ctrl
  )

  X <- model.matrix(
    ~ x1 + x2,
    data = dat
  )

  Z <- model.matrix(
    ~ z1 + z2 - 1,
    data = dat
  )

  fit_internal <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = X,
    Z = Z,
    family = binomial(),
    M = 4,
    lambda = 1,
    maxit = 30,
    tol = 1e-5,
    trace = FALSE
  )

  expect_equal(
    fit_public$alpha_hat,
    fit_internal$alpha_hat
  )

  expect_equal(
    fit_public$beta_hat,
    fit_internal$beta_hat
  )

  expect_equal(
    fit_public$eta_hat,
    fit_internal$eta_hat
  )

  expect_equal(
    fit_public$fitted_mean,
    fit_internal$fitted_mean
  )

  expect_equal(
    fit_public$logLik,
    fit_internal$logLik
  )

  expect_equal(
    fit_public$converged,
    fit_internal$converged
  )
})


test_that("formula and matrix interfaces give equivalent Poisson fits", {

  set.seed(2027)

  n <- 200

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n),
    exposure = runif(n, 100, 1000)
  )

  eta <- -5 +
    0.4 * dat$x1 -
    0.2 * dat$x2 +
    0.4 * sin(dat$z1 + dat$z2)

  dat$y <- rpois(
    n,
    lambda = dat$exposure * exp(eta)
  )

  ctrl <- pgplsim_control(
    maxit = 30,
    tol = 1e-5,
    trace = FALSE
  )

  fit_public <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = log(dat$exposure),
    M = 4,
    lambda = 1,
    control = ctrl
  )

  X <- model.matrix(
    ~ x1 + x2,
    data = dat
  )

  Z <- model.matrix(
    ~ z1 + z2 - 1,
    data = dat
  )

  fit_internal <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = X,
    Z = Z,
    family = poisson(),

    # Internal engine uses positive exposure scale.
    offset = dat$exposure,

    M = 4,
    lambda = 1,
    maxit = 30,
    tol = 1e-5,
    trace = FALSE
  )

  expect_equal(
    fit_public$alpha_hat,
    fit_internal$alpha_hat
  )

  expect_equal(
    fit_public$beta_hat,
    fit_internal$beta_hat
  )

  expect_equal(
    fit_public$eta_hat,
    fit_internal$eta_hat
  )

  expect_equal(
    fit_public$fitted_mean,
    fit_internal$fitted_mean
  )

  expect_equal(
    fit_public$logLik,
    fit_internal$logLik
  )

  expect_equal(
    fit_public$converged,
    fit_internal$converged
  )

  # Explicitly lock down the public/internal offset distinction.
  expect_equal(
    fit_public$offset,
    log(dat$exposure)
  )

  expect_equal(
    fit_public$exposure,
    dat$exposure
  )
})


test_that("formula frontend constructs the same X and Z matrices", {

  set.seed(2028)

  n <- 120

  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  X_expected <- model.matrix(
    ~ x1 + x2,
    data = dat
  )

  Z_expected <- model.matrix(
    ~ z1 + z2 - 1,
    data = dat
  )

  expect_equal(
    fit$X,
    X_expected
  )

  expect_equal(
    fit$Z,
    Z_expected
  )

  expect_equal(
    fit$y,
    dat$y
  )
})
