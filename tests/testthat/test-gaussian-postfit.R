# Postfit Gaussian validation against the July 30, 2026 Lu method.
make_step16_gaussian_data <- function() {
  set.seed(20261123)
  n <- 240L
  d <- data.frame(x1 = rnorm(n), z1 = runif(n, -1, 1),
                  z2 = runif(n, -1, 1), z3 = runif(n, -1, 1))
  u <- as.vector(as.matrix(d[c("z1", "z2", "z3")]) %*%
                   rep(1 / sqrt(3), 3))
  f <- sin(2 * u)
  d$y <- 0.8 * d$x1 + f - mean(f) + rnorm(n, sd = 0.4)
  d
}

make_step16_gaussian_fit <- function() {
  d <- make_step16_gaussian_data()
  fit <- pgplsim(y ~ x1, index = ~ z1 + z2 + z3,
                 data = d, family = gaussian(), M = 5,
                 control = pgplsim_control(maxit = 60, tol = 1e-5))
  list(fit = fit, data = d)
}

test_that("Gaussian prediction is identity-linked and decomposes correctly", {
  z <- make_step16_gaussian_fit()
  f <- z$fit
  d <- z$data
  expect_equal(as.numeric(predict(f, type = "response")),
               as.numeric(f$fitted_mean), tolerance = 1e-6)
  expect_equal(as.numeric(predict(f, type = "link")),
               as.numeric(f$fitted_mean), tolerance = 1e-6)
  expect_equal(as.numeric(predict(f, type = "response")),
               as.numeric(predict(f, type = "linear") +
                          predict(f, type = "smooth")), tolerance = 1e-6)
  newdata <- d[1:15, , drop = FALSE]
  expect_equal(as.numeric(predict(f, newdata = newdata)),
               as.numeric(predict(f, newX = f$X[1:15, , drop = FALSE],
                                  newZ = f$Z[1:15, , drop = FALSE])),
               tolerance = 1e-6)
  expect_error(predict(f, offset = rep(0, nrow(d))), "Poisson")
})

test_that("Gaussian summary reports adjusted dispersion and normal reference", {
  z <- make_step16_gaussian_fit()
  f <- z$fit
  s <- summary(f)
  expect_s3_class(s, "summary.pgplsim")
  expect_equal(s$dispersion, f$dispersion)
  expect_equal(s$residual_df, f$residual_df)
  expect_equal(s$effective_df, f$effective_df)
  expect_equal(s$rss, f$rss)
  expect_equal(s$logLik_mle, f$logLik_mle)
  expect_true(all(c("Estimate", "SE", "z", "p") %in% names(s$linear)))
  expect_equal(s$linear$z, s$linear$Estimate / s$linear$SE)
  expect_equal(s$linear$p,
               2 * stats::pnorm(abs(s$linear$z), lower.tail = FALSE))
  expect_output(print(s), "Residual dispersion")
})

test_that("Gaussian case bootstrap refits and yields finite SEs", {
  z <- make_step16_gaussian_fit()
  f <- z$fit
  args <- list(object = f, B = 3L, method = "case", seed = 2026,
               max_attempts = 30L, trace = FALSE)
  boot1 <- do.call(bootstrap_pgplsim, args)
  boot2 <- do.call(bootstrap_pgplsim, args)
  expect_s3_class(boot1, "bootstrap.pgplsim")
  expect_equal(boot1$B, 3L)
  expect_equal(nrow(boot1$estimates_alpha), 3)
  expect_equal(nrow(boot1$estimates_beta), 3)
  expect_true(all(is.finite(boot1$se_alpha)))
  expect_true(all(is.finite(boot1$se_beta)))
  expect_equal(boot1$estimates_alpha, boot2$estimates_alpha)
  expect_equal(boot1$estimates_beta, boot2$estimates_beta)
  expect_equal(vcov(f, type = "bootstrap", bootstrap = boot1),
               stats::cov(boot1$estimates_parameters))
  expect_equal(boot1$failures, boot1$attempts - boot1$B)
  expect_error(bootstrap_pgplsim(f, B = 2, method = "stratified"),
               "method = 'case'")
  expect_error(bootstrap_pgplsim(f, B = 2, offset = rep(0, nrow(z$data))),
               "Gaussian bootstrap")
})
