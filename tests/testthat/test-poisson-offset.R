test_that("Poisson likelihood and offsets are consistent", {
  set.seed(20261106)
  n <- 180
  dat <- data.frame(
    x1 = rnorm(n),
    z1 = runif(n, -1, 1),
    z2 = runif(n, -1, 1),
    exposure = runif(n, 0.5, 5)
  )
  u <- 0.8 * dat$z1 + 0.6 * dat$z2
  eta_true <- 0.2 + 0.5 * dat$x1 + 0.3 * sin(2 * u)
  dat$y <- rpois(n, dat$exposure * exp(eta_true))
  off <- log(dat$exposure)

  fit <- pgplsim(
    y ~ x1, index = ~ z1 + z2, data = dat,
    family = poisson(), offset = off, M = 4,
    control = pgplsim_control(maxit = 100, tol = 1e-6)
  )

  mu <- as.numeric(fit$fitted_mean)
  expect_equal(
    mu, exp(as.numeric(fit$fitted_link) + off),
    tolerance = 1e-7
  )
  expect_equal(
    as.numeric(fit$logLik),
    sum(dat$y * log(mu) - mu),
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(predict(fit, type = "response")),
    mu, tolerance = 1e-7
  )

  mu_new <- as.numeric(predict(
    fit, newdata = dat, type = "response", offset = off
  ))
  expect_equal(mu_new, mu, tolerance = 1e-6)

  mu_double <- as.numeric(predict(
    fit, newdata = dat, type = "response",
    offset = off + log(2)
  ))
  expect_equal(mu_double, 2 * mu_new, tolerance = 1e-6)

  expect_equal(
    as.numeric(predict(fit, type = "link")),
    log(mu), tolerance = 1e-7
  )
})
