test_that("estimated alpha satisfies identifiability constraints", {

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

  expect_equal(
    sum(fit$alpha_hat^2),
    1,
    tolerance = 1e-8
  )

  expect_gte(
    fit$alpha_hat[1],
    0
  )
})
