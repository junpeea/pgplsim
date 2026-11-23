# Gaussian covariance validation, Lu (2026), Sections 4--5.
# This file reconstructs the four covariance formulas independently.

make_gaussian_covariance_lu_fit <- function() {
  set.seed(20261123)
  n <- 220L
  dat <- data.frame(x1 = rnorm(n), z1 = runif(n, -1, 1),
                    z2 = runif(n, -1, 1), z3 = runif(n, -1, 1))
  u <- (dat$z1 + dat$z2 + dat$z3) / sqrt(3)
  dat$y <- 0.6 + 0.8 * dat$x1 + sin(2 * u) + rnorm(n, sd = 0.5)
  pgplsim(y ~ x1, index = ~ z1 + z2 + z3,
          data = dat, family = gaussian(), M = 5,
          control = pgplsim_control(maxit = 45, tol = 1e-6))
}

test_that("Gaussian model and sandwich covariance match direct matrix formulas", {
  fit <- make_gaussian_covariance_lu_fit()
  E <- as.matrix(fit$estimating_design)
  d <- nrow(E)
  P <- pgplsim:::generate_P_matrix_pkg(d, length(fit$gamma_hat), 2)
  Kinv <- solve(tcrossprod(E) + fit$lambda * P)
  V_model <- fit$dispersion * Kinv
  score <- as.matrix(fit$score_contributions)
  H <- crossprod(score)
  V_sandwich <- V_model %*% H %*% V_model
  expect_equal(unname(vcov(fit, type = "model", component = "full")),
               unname(V_model), tolerance = 1e-5)
  expect_equal(unname(vcov(fit, type = "sandwich", component = "full")),
               unname(V_sandwich), tolerance = 1e-5)
})

test_that("Gaussian leverage correction uses squared inverse 1-h", {
  fit <- make_gaussian_covariance_lu_fit()
  E <- as.matrix(fit$estimating_design)
  score <- as.matrix(fit$score_contributions)
  h <- colSums(E * solve(tcrossprod(E), E))
  expect_true(all(h >= -1e-8 & h < 1))
  expect_gt(max(h), 0)
  bread <- unname(vcov(fit, type = "model", component = "full"))
  adjusted <- score / (1 - h)
  V_direct <- bread %*% crossprod(adjusted) %*% bread
  V_pkg <- unname(vcov(fit, type = "leverage", component = "full"))
  expect_equal(V_pkg, V_direct, tolerance = 1e-5)
  # It must not accidentally reduce to the legacy HC2 adjustment.
  V_hc2 <- bread %*% crossprod(score / sqrt(1 - h)) %*% bread
  expect_gt(max(abs(V_pkg - V_hc2)), 1e-9)
})

test_that("Gaussian finite correction uses the Morel additive trace term", {
  fit <- make_gaussian_covariance_lu_fit()
  bread <- unname(vcov(fit, type = "model", component = "full"))
  H <- crossprod(as.matrix(fit$score_contributions))
  V_robust <- bread %*% H %*% bread
  n <- length(fit$y)
  rho <- fit$effective_df
  delta <- min(0.5, rho / (n - rho))
  kappa <- max(1, sum(diag(bread %*% H)) / rho)
  V_direct <- V_robust + delta * kappa * bread
  V_pkg <- unname(vcov(fit, type = "finite", component = "full"))
  expect_equal(V_pkg, V_direct, tolerance = 1e-5)
  expect_true(all(eigen(V_pkg, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(unname(vcov(fit, type = "finite", component = "full",
                           finite_df = rho)), V_direct, tolerance = 1e-5)
  expect_error(vcov(fit, type = "finite", finite_df = n), "between 0 and n")
})

test_that("Gaussian eta-to-alpha delta transform retains cross covariance", {
  fit <- make_gaussian_covariance_lu_fit()
  eta <- as.numeric(fit$eta_hat)
  h <- 1 + sum(eta^2)
  J <- rbind(-eta / h^(3/2),
             diag(length(eta)) / sqrt(h) - outer(eta, eta) / h^(3/2))
  p_eta <- length(eta)
  p_beta <- length(fit$beta_hat)
  d <- nrow(as.matrix(fit$estimating_design))
  T <- matrix(0, nrow = p_eta + 1L + p_beta, ncol = d)
  T[seq_len(p_eta + 1L), seq_len(p_eta)] <- J
  T[p_eta + 1L + seq_len(p_beta), p_eta + seq_len(p_beta)] <- diag(p_beta)
  for (kind in c("model", "sandwich", "leverage", "finite")) {
    full <- unname(vcov(fit, type = kind, component = "full"))
    expected <- T %*% full %*% t(T)
    observed <- unname(vcov(fit, type = kind, component = "parameters"))
    expect_equal(observed, expected, tolerance = 1e-6)
    expect_equal(unname(vcov(fit, type = kind, component = "alpha")),
                 expected[seq_len(p_eta + 1L), seq_len(p_eta + 1L)],
                 tolerance = 1e-6)
  }
})
