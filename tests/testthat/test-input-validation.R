test_that("incompatible dimensions produce an error", {

  y <- rbinom(100, 1, 0.5)
  X <- matrix(rnorm(200), nrow = 100)
  Z <- matrix(rnorm(297), nrow = 99)

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y,
      X = X,
      Z = Z,
      family = binomial()
    ),
    "incompatible"
  )
})

test_that("binomial outcomes must be zero or one", {

  y <- c(0, 1, 2, 0)
  X <- matrix(rnorm(8), nrow = 4)
  Z <- matrix(rnorm(12), nrow = 4)

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y,
      X = X,
      Z = Z,
      family = binomial()
    ),
    "must contain only 0 and 1"
  )
})

test_that("Poisson exposure must be positive", {

  dat <- make_poisson_test_data(n = 50)

  bad_exposure <- dat$exposure
  bad_exposure[1] <- 0

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = bad_exposure
    ),
    "must be positive"
  )
})


test_that("missing and non-finite responses are rejected", {

  dat <- make_binomial_test_data(n = 60, seed = 1201)

  y_na <- dat$y
  y_na[1] <- NA_real_

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y_na,
      X = dat$X,
      Z = dat$Z,
      family = binomial()
    ),
    "finite, non-missing"
  )

  y_inf <- as.numeric(dat$y)
  y_inf[1] <- Inf

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y_inf,
      X = dat$X,
      Z = dat$Z,
      family = binomial()
    ),
    "finite, non-missing"
  )

  y_nan <- as.numeric(dat$y)
  y_nan[1] <- NaN

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y_nan,
      X = dat$X,
      Z = dat$Z,
      family = binomial()
    ),
    "finite, non-missing"
  )
})


test_that("missing and non-finite design values are rejected", {

  dat <- make_binomial_test_data(n = 60, seed = 1202)

  X_na <- dat$X
  X_na[1, 1] <- NA_real_

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = X_na,
      Z = dat$Z,
      family = binomial()
    )
  )

  X_inf <- dat$X
  X_inf[1, 1] <- Inf

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = X_inf,
      Z = dat$Z,
      family = binomial()
    )
  )

  Z_na <- dat$Z
  Z_na[1, 1] <- NA_real_

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = Z_na,
      family = binomial()
    )
  )

  Z_nan <- dat$Z
  Z_nan[1, 1] <- NaN

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = Z_nan,
      family = binomial()
    )
  )
})


test_that("Poisson responses cannot be negative or non-integer", {

  dat <- make_poisson_test_data(n = 60, seed = 1203)

  y_negative <- dat$y
  y_negative[1] <- -1

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y_negative,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = dat$exposure
    ),
    "non-negative integer counts"
  )

  y_fractional <- as.numeric(dat$y)
  y_fractional[1] <- y_fractional[1] + 0.5

  expect_error(
    pgplsim:::pgplsim_fit(
      y = y_fractional,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = dat$exposure
    ),
    "non-negative integer counts"
  )
})


test_that("Poisson data may contain zero counts", {

  dat <- make_poisson_test_data(n = 100, seed = 1204)
  dat$y[seq_len(10)] <- 0

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = poisson(),
    offset = dat$exposure,
    M = 4,
    maxit = 30
  )

  expect_s3_class(fit, "pgplsim")
  expect_true(all(fit$fitted_mean > 0))
})


test_that("offset must have the correct length", {

  dat <- make_poisson_test_data(n = 60, seed = 1205)

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = dat$exposure[-1]
    ),
    "same length as y"
  )
})


test_that("missing and non-finite Poisson offsets are rejected", {

  dat <- make_poisson_test_data(n = 60, seed = 1206)

  offset_na <- dat$exposure
  offset_na[1] <- NA_real_

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = offset_na
    )
  )

  offset_inf <- dat$exposure
  offset_inf[1] <- Inf

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = offset_inf
    )
  )
})


test_that("unsupported families are rejected", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 123
  )

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = gaussian()
    ),
    "Only `binomial()` and `poisson()` are currently supported",
    fixed = TRUE
  )
})


test_that("constant X columns are handled without crashing", {

  set.seed(2026)

  n <- 80
  X <- cbind(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  Z <- cbind(
    Z1 = rnorm(n),
    Z2 = rnorm(n),
    Z3 = rnorm(n)
  )

  eta <- 0.5 * X[, 1] - 0.3 * X[, 2] +
    0.6 * Z[, 1]

  y <- rbinom(
    n,
    size = 1,
    prob = plogis(eta)
  )

  X[, 1] <- 1

  fit <- pgplsim:::pgplsim_fit(
    y = y,
    X = X,
    Z = Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  expect_s3_class(fit, "pgplsim")
  expect_true(is.logical(fit$converged))
})


test_that("constant Z columns are rejected or fail safely", {

  dat <- make_binomial_test_data(n = 80, seed = 1209)
  Z_constant <- cbind(dat$Z, constant = 1)

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = Z_constant,
      family = binomial(),
      M = 4,
      maxit = 20
    )
  )
})


test_that("rank-deficient X is rejected or fails safely", {

  dat <- make_binomial_test_data(n = 80, seed = 1210)
  X_rank_deficient <- cbind(
    dat$X,
    duplicate = dat$X[, 1]
  )

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = X_rank_deficient,
      Z = dat$Z,
      family = binomial(),
      M = 4,
      maxit = 20
    )
  )
})


test_that("very small samples fail safely", {

  dat <- make_binomial_test_data(n = 5, seed = 1211)

  expect_error(
    pgplsim:::pgplsim_fit(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = binomial(),
      M = 4,
      maxit = 10
    )
  )
})


test_that("all-zero Poisson responses are handled safely", {

  set.seed(2026)

  n <- 80
  X <- cbind(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  Z <- cbind(
    Z1 = rnorm(n),
    Z2 = rnorm(n),
    Z3 = rnorm(n)
  )

  exposure <- rep(1, n)
  y <- rep(0, n)

  fit <- pgplsim:::pgplsim_fit(
    y = y,
    X = X,
    Z = Z,
    family = poisson(),
    offset = exposure,
    M = 4,
    maxit = 30
  )

  expect_s3_class(fit, "pgplsim")
  expect_true(is.logical(fit$converged))
})
