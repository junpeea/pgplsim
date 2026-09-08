test_that("public formula interface fits binomial data", {

  set.seed(2026)

  n <- 150

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  eta <- 0.5 * dat$x1 -
    0.3 * dat$x2 +
    sin(dat$z1 + dat$z2)

  dat$y <- rbinom(
    n,
    size = 1,
    prob = plogis(eta)
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4,
    control = pgplsim_control(
      maxit = 30,
      tol = 1e-5
    )
  )

  expect_s3_class(fit, "pgplsim")
  expect_true(fit$converged)

  expect_equal(
    names(fit$beta_hat),
    c("(Intercept)", "x1", "x2")
  )

  expect_equal(
    names(fit$alpha_hat),
    c("z1", "z2")
  )

  expect_equal(nrow(fit$X), n)
  expect_equal(nrow(fit$Z), n)
  expect_equal(fit$nobs, n)
})


test_that("formula and index metadata are stored", {

  set.seed(2026)

  n <- 100

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

  expect_equal(
    fit$formula,
    y ~ x1 + x2
  )

  expect_equal(
    fit$index,
    ~ z1 + z2
  )

  expect_true(inherits(fit$terms, "terms"))
  expect_true(inherits(fit$index_terms, "terms"))

  expect_equal(
    colnames(fit$X),
    c("(Intercept)", "x1", "x2")
  )

  expect_equal(
    colnames(fit$Z),
    c("z1", "z2")
  )

  expect_false(is.null(fit$call))
})


test_that("public Poisson offset uses link-scale semantics", {

  set.seed(2026)

  n <- 200

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n),
    exposure = runif(n, 100, 1000)
  )

  eta <- -5 +
    0.3 * dat$x1 -
    0.2 * dat$x2 +
    0.4 * sin(dat$z1 + dat$z2)

  mu <- dat$exposure * exp(eta)

  dat$y <- rpois(n, mu)

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = log(dat$exposure),
    M = 4,
    control = pgplsim_control(
      maxit = 30
    )
  )

  expect_s3_class(fit, "pgplsim")

  expect_true(fit$offset_supplied)

  expect_equal(
    fit$offset,
    log(dat$exposure)
  )

  expect_equal(
    fit$exposure,
    dat$exposure
  )
})


test_that("public Poisson interface agrees with internal exposure interface", {

  set.seed(2026)

  n <- 150

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
    0.3 * sin(dat$z1 + dat$z2)

  dat$y <- rpois(
    n,
    dat$exposure * exp(eta)
  )

  fit_public <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = log(dat$exposure),
    M = 4,
    control = pgplsim_control(
      maxit = 30
    )
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
    offset = dat$exposure,
    M = 4,
    maxit = 30
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
    fit_public$fitted_mean,
    fit_internal$fitted_mean
  )

  expect_equal(
    fit_public$logLik,
    fit_internal$logLik
  )
})


test_that("missing values use one synchronized analysis sample", {

  set.seed(2026)

  n <- 100

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n),
    exposure = runif(n, 100, 1000)
  )

  dat$y <- rpois(n, 2)

  dat$x1[2] <- NA
  dat$z1[5] <- NA

  off <- log(dat$exposure)
  off[8] <- NA

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = off,
    M = 4
  )

  expect_equal(fit$nobs, 97)

  expect_equal(nrow(fit$X), 97)
  expect_equal(nrow(fit$Z), 97)
  expect_equal(length(fit$y), 97)
  expect_equal(length(fit$offset), 97)
  expect_equal(length(fit$exposure), 97)

  expect_equal(
    as.integer(fit$na.action),
    c(2L, 5L, 8L)
  )
})


test_that("factor predictors are handled by the linear formula", {

  set.seed(2026)

  n <- 150

  dat <- data.frame(
    x1 = rnorm(n),
    group = factor(
      sample(c("A", "B", "C"), n, replace = TRUE)
    ),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  dat$y <- rbinom(n, 1, 0.5)

  fit <- pgplsim(
    y ~ x1 + group,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  expect_s3_class(fit, "pgplsim")

  expect_true("(Intercept)" %in% colnames(fit$X))
  expect_true("x1" %in% colnames(fit$X))

  expect_equal(
    ncol(fit$X),
    4
  )

  expect_true(
    all(c("z1", "z2") %in% colnames(fit$Z))
  )
})


test_that("nonnumeric index variables are rejected", {

  set.seed(2026)

  n <- 100

  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    z1 = rnorm(n),
    group = factor(
      sample(c("A", "B"), n, replace = TRUE)
    )
  )

  expect_error(
    pgplsim(
      y ~ x1,
      index = ~ z1 + group,
      data = dat,
      family = binomial()
    ),
    "numeric"
  )
})


test_that("binomial models reject offsets", {

  set.seed(2026)

  n <- 100

  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  expect_error(
    pgplsim(
      y ~ x1,
      index = ~ z1 + z2,
      data = dat,
      family = binomial(),
      offset = rep(0, n)
    ),
    "only for Poisson"
  )
})


test_that("unsupported families are rejected by public interface", {

  set.seed(2026)

  n <- 100

  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  expect_error(
    pgplsim(
      y ~ x1,
      index = ~ z1 + z2,
      data = dat,
      family = gaussian()
    ),
    "binomial"
  )
})


test_that("pgplsim_control settings reach fitted object", {

  set.seed(2026)

  n <- 120

  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  ctrl <- pgplsim_control(
    maxit = 35,
    tol = 1e-6,
    trace = FALSE
  )

  fit <- pgplsim(
    y ~ x1,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4,
    control = ctrl
  )

  expect_equal(fit$control$maxit, 35)
  expect_equal(fit$control$tol, 1e-6)
  expect_false(fit$control$trace)
})
