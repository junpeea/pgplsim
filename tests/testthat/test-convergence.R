test_that("convergence flag is TRUE when tolerance is easily satisfied", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1101
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = poisson(),
    offset = dat$exposure,
    M = 4,
    maxit = 30,
    tol = 1e-5
  )

  expect_true(fit$converged)
  expect_lte(fit$iterations, 5)
  expect_true(is.finite(fit$convergence_difference))
})


test_that("non-convergence is recorded when iteration limit is reached", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1102
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 1,
    tol = .Machine$double.eps
  )

  expect_false(fit$converged)
  expect_equal(fit$iterations, 1)
  expect_true(is.finite(fit$convergence_difference))
})
