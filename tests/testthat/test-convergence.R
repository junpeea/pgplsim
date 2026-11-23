test_that("convergence is reached within maxit at the requested tolerance", {

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
  expect_gte(fit$iterations, 1L)
  expect_lt(fit$iterations, 30L)
  expect_true(is.finite(fit$convergence_difference))
  expect_lt(fit$convergence_difference, 1e-5)
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
