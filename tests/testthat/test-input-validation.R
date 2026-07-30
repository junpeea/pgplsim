test_that("incompatible dimensions produce an error", {

  y <- rbinom(100, 1, 0.5)
  X <- matrix(rnorm(200), nrow = 100)
  Z <- matrix(rnorm(297), nrow = 99)

  expect_error(
    LuGPLSIM(
      y = y,
      X = X,
      Z = Z,
      family = binomial()
    ),
    "Dimensions of y, X, and Z are incompatible"
  )
})

test_that("binomial outcomes must be zero or one", {

  y <- c(0, 1, 2, 0)
  X <- matrix(rnorm(8), nrow = 4)
  Z <- matrix(rnorm(12), nrow = 4)

  expect_error(
    LuGPLSIM(
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
    LuGPLSIM(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      family = poisson(),
      offset = bad_exposure
    ),
    "must be positive"
  )
})
