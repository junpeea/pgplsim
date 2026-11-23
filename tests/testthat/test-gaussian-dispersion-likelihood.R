# Independent Gaussian diagnostics for Phase 1, Step 1.3.
# These checks are intentionally based on the returned means, scores and
# Fisher-design matrix rather than reusing findLogLik_gaussian_pkg().

test_that("Gaussian variance uses residual, not model, degrees of freedom", {
  set.seed(20261124)
  n <- 220L
  dat <- data.frame(x1 = rnorm(n), z1 = runif(n, -1, 1),
                    z2 = runif(n, -1, 1), z3 = runif(n, -1, 1))
  u <- (dat$z1 + dat$z2 + dat$z3) / sqrt(3)
  dat$y <- 1 + 0.7 * dat$x1 + sin(2 * u) + rnorm(n, sd = 0.6)

  fit <- pgplsim(y ~ x1, index = ~ z1 + z2 + z3,
                 data = dat, family = gaussian(), M = 5,
                 control = pgplsim_control(maxit = 45, tol = 1e-6))

  rss <- sum((dat$y - as.numeric(fit$fitted_mean))^2)
  expect_equal(fit$rss, rss, tolerance = 1e-7)
  expect_equal(fit$effective_df + fit$residual_df, n, tolerance = 1e-8)
  expect_gt(fit$residual_df, 0)
  expect_gt(fit$effective_df, 0)
  expect_lt(fit$effective_df, n)
  expect_equal(fit$dispersion, rss / fit$residual_df, tolerance = 1e-8)
  expect_equal(fit$dispersion_mle, rss / n, tolerance = 1e-8)

  # Independently reconstruct the effective MODEL df from unit-dispersion
  # estimating derivatives and the exact penalty multiplier stored in fit.
  E <- as.matrix(fit$estimating_design)  # d x n, internally index sorted
  d <- nrow(E)
  P <- pgplsim:::generate_P_matrix_pkg(d, length(fit$gamma_hat), 2)
  K <- tcrossprod(E) + fit$lambda * P
  edf_direct <- sum(diag(solve(K, tcrossprod(E))))
  expect_equal(fit$effective_df, edf_direct, tolerance = 1e-5)

  # Check the Gaussian covariance on the full unconstrained parameter scale.
  # The manuscript's Fisher information includes 1/sigma^2; hence V = sigma^2 K^-1.
  expect_equal(unname(fit$vcov_full), unname(fit$dispersion * solve(K)),
               tolerance = 1e-5)
})

test_that("Gaussian full and penalized likelihoods have explicit scaling", {
  set.seed(20261125)
  n <- 180L
  dat <- data.frame(x = rnorm(n), z1 = runif(n), z2 = runif(n))
  dat$y <- 0.4 + 0.8 * dat$x + sin(2 * (dat$z1 + dat$z2)) +
    rnorm(n, sd = 0.45)
  fit <- pgplsim(y ~ x, index = ~ z1 + z2, data = dat,
                 family = gaussian(), M = 5)

  mu <- as.numeric(fit$fitted_mean)
  rss <- sum((dat$y - mu)^2)
  lik_pearson <- sum(dnorm(dat$y, mean = mu,
                            sd = sqrt(rss / fit$residual_df), log = TRUE))
  lik_mle <- sum(dnorm(dat$y, mean = mu,
                        sd = sqrt(rss / n), log = TRUE))
  expect_equal(as.numeric(fit$logLik), lik_pearson, tolerance = 1e-7)
  expect_equal(as.numeric(fit$logLik_mle), lik_mle, tolerance = 1e-7)
  expect_gte(fit$logLik_mle + 1e-7, fit$logLik)

  # The manuscript uses lambda_manuscript^2. The package stores its multiplier.
  DTD <- pgplsim:::generate_P_matrix_pkg(
    length(fit$eta_hat) + length(fit$beta_hat) + length(fit$gamma_hat),
    length(fit$gamma_hat), 2)
  theta <- c(unname(fit$eta_hat), unname(fit$beta_hat),
             unname(fit$gamma_hat))
  roughness <- as.numeric(crossprod(theta, DTD %*% theta))
  penalty <- 0.5 * fit$lambda * roughness
  expect_equal(fit$penalty_quadratic, roughness, tolerance = 1e-7)
  expect_equal(fit$penalty_value, penalty, tolerance = 1e-7)
  expect_equal(fit$penalized_logLik,
               fit$logLik - penalty / fit$dispersion, tolerance = 1e-7)

  # Check the observation-level Gaussian score scaling independently.
  # The likelihood score per observation is (y - mu) E_i / sigma^2.
  # Compare in original input order only up to permutation: E is index sorted.
  E <- t(as.matrix(fit$estimating_design))
  score <- as.matrix(fit$score_contributions)
  expect_equal(dim(score), dim(E))
  expect_equal(unname(crossprod(score)),
               unname(crossprod(E * as.vector(
                 fit$residuals[order(as.vector(fit$Z %*% fit$alpha_hat))]) /
                 fit$dispersion)), tolerance = 1e-5)
})
