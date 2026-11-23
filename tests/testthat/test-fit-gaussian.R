test_that("Gaussian identity-link fit returns finite structural results", {
  set.seed(20261123)
  n <- 240
  dat <- data.frame(
    x1 = rnorm(n),
    z1 = runif(n, -1, 1),
    z2 = runif(n, -1, 1),
    z3 = runif(n, -1, 1)
  )
  alpha0 <- rep(1 / sqrt(3), 3)
  u <- as.vector(as.matrix(dat[c("z1", "z2", "z3")]) %*% alpha0)
  dat$y <- 0.4 + 0.8 * dat$x1 + sin(2 * u) + rnorm(n, sd = 0.35)

  fit <- pgplsim(
    y ~ x1, index = ~ z1 + z2 + z3,
    data = dat, family = gaussian(), M = 5,
    control = pgplsim_control(maxit = 40, tol = 1e-5)
  )

  expect_s3_class(fit, "pgplsim")
  expect_equal(length(fit$fitted_mean), n)
  expect_equal(sum(fit$alpha_hat^2), 1, tolerance = 1e-8)
  expect_true(fit$alpha_hat[1] > 0)
  expect_true(is.finite(fit$dispersion) && fit$dispersion > 0)
  expect_true(is.finite(fit$logLik))
  expect_true(is.finite(fit$lambda) && fit$lambda > 0)
  expect_equal(fit$fitted_link, fit$fitted_mean, tolerance = 1e-9)
  expect_equal(mean(fit$smooth_component), 0, tolerance = 1e-8)
  expect_equal(
    as.numeric(predict(fit, type = "response")),
    as.numeric(fit$fitted_mean), tolerance = 1e-6
  )
  expect_equal(
    as.numeric(predict(fit, type = "link")),
    as.numeric(predict(fit, type = "response")), tolerance = 1e-9
  )
  expect_equal(dim(vcov(fit, type = "model", component = "beta")),
               c(ncol(fit$X), ncol(fit$X)))
})

test_that("Gaussian nonidentity links and offsets are explicitly rejected", {
  set.seed(13)
  d <- data.frame(y = rnorm(60), x = rnorm(60),
                  z1 = runif(60), z2 = runif(60))
  expect_error(pgplsim(y ~ x, index = ~ z1 + z2, data = d,
                       family = gaussian(link = "log")), "identity")
  expect_error(pgplsim(y ~ x, index = ~ z1 + z2, data = d,
                       family = gaussian(), offset = rep(0, 60)), "Poisson")
})
