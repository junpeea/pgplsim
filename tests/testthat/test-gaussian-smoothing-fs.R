test_that("Second-order P-spline penalty has the correct null space", {
  q <- 6L
  P <- pgplsim:::generate_P_matrix_pkg(q + 3L, q, 2)
  Pg <- P[4:9, 4:9, drop = FALSE]
  D <- matrix(0, nrow = q - 2L, ncol = q)
  for (j in seq_len(q - 2L)) D[j, j:(j + 2L)] <- c(1, -2, 1)
  expect_equal(Pg, crossprod(D), tolerance = 1e-12)
  expect_equal(as.numeric(crossprod(rep(1, q), Pg %*% rep(1, q))),
               0, tolerance = 1e-12)
  expect_equal(as.numeric(crossprod(seq_len(q), Pg %*% seq_len(q))),
               0, tolerance = 1e-12)
  expect_gt(as.numeric(crossprod(seq_len(q)^2,
                                  Pg %*% seq_len(q)^2)), 0)
})

# Phase 1 Step 1.4: independent matrix checks of Lu (2026), Equation (9).

test_that("Gaussian FS multiplier matches manuscript pseudoinverse formula", {
  set.seed(14)
  npar <- 9L
  ngamma <- 5L
  A <- matrix(rnorm(npar * npar), npar, npar)
  s <- 1.7  # package multiplier, i.e. manuscript lambda^2
  P <- pgplsim:::generate_P_matrix_pkg(npar, ngamma, 2)
  K <- crossprod(A) + diag(2, npar) + s * P
  gamma <- c(0.6, -0.2, 0.7, -0.1, 0.8)
  mfit <- list(H.inv = -solve(K),
               estimating_design = matrix(0, npar, 20))
  actual <- pgplsim:::update_lambda_gaussian_FS_pkg(mfit, gamma, s)
  # Direct Eq. (9) with P_s = s P and generalized inverse of P_s.
  numerator <- sum(diag(MASS::ginv(s * P) %*% P)) -
    sum(diag(solve(K) %*% P))
  theta <- c(rep(0, npar - ngamma), gamma)
  denominator <- as.numeric(crossprod(theta, P %*% theta))
  expected <- s * numerator / denominator
  expect_true(actual$updated)
  expect_equal(actual$lambda, expected, tolerance = 1e-9)
  expect_equal(actual$penalty_rank, ngamma - 2L)
  expect_gt(actual$lambda, 0)
})

test_that("Gaussian FS retains lambda for zero roughness", {
  P <- pgplsim:::generate_P_matrix_pkg(8, 4, 2)
  mfit <- list(H.inv = -solve(diag(8) + P),
               estimating_design = matrix(0, 8, 10))
  ans <- pgplsim:::update_lambda_gaussian_FS_pkg(
    mfit, gamma = rep(1, 4), lambda = 1)
  expect_false(ans$updated)
  expect_equal(ans$lambda, 1)
  expect_error(pgplsim:::update_lambda_gaussian_FS_pkg(
    mfit, gamma = rep(1, 4), lambda = 0), "positive")
})

test_that("Gaussian smoothing uses fixed-lambda inner and outer updates", {
  set.seed(20261123)
  n <- 160L
  d <- data.frame(x1 = rnorm(n), z1 = runif(n, -1, 1),
                  z2 = runif(n, -1, 1), z3 = runif(n, -1, 1))
  u <- (d$z1 + d$z2 + d$z3) / sqrt(3)
  d$y <- 0.3 + 0.7 * d$x1 + sin(2 * u) + rnorm(n, sd = 0.45)
  fit <- pgplsim(y ~ x1, index = ~ z1 + z2 + z3,
                 data = d, family = gaussian(), M = 5,
                 control = pgplsim_control(maxit = 40, tol = 1e-5))
  sm <- fit$smoothing_diagnostics
  expect_true(is.list(sm))
  expect_equal(length(sm$lambda_history), fit$outer_iterations)
  expect_equal(length(sm$proposed_lambda_history), fit$outer_iterations)
  expect_equal(length(sm$inner_iterations), fit$outer_iterations)
  expect_gt(min(sm$lambda_history), 0)
  expect_true(all(is.finite(sm$proposed_lambda_history)))
  expect_equal(fit$lambda, tail(sm$lambda_history, 1), tolerance = 1e-12)
  if (fit$outer_iterations > 1L) {
    expect_equal(sm$lambda_history[-1L],
                 head(sm$proposed_lambda_history, -1L), tolerance = 1e-12)
  }
  # At the final fitted state, Eq. (9) is reproducible without refitting.
  endfit <- pgplsim:::findLogLik_gaussian_pkg(
    Y = fit$y, X = fit$X, Z = fit$Z,
    alpha = unname(fit$alpha_hat), beta = unname(fit$beta_hat),
    gamma = unname(fit$gamma_hat), lambda = fit$lambda, qn = fit$qn)
  upd <- pgplsim:::update_lambda_gaussian_FS_pkg(
    endfit, unname(fit$gamma_hat), fit$lambda)
  expect_equal(tail(sm$proposed_lambda_history, 1),
               upd$lambda, tolerance = 1e-7)
})
