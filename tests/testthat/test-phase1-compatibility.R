# Phase 1 final three-family compatibility regression tests.
# Designed to augment, not replace, the individual family tests.

make_phase1_compatibility_data <- function(n = 200L) {
  set.seed(20261123)
  d <- data.frame(
    x1 = rnorm(n),
    z1 = runif(n, -1, 1),
    z2 = runif(n, -1, 1),
    z3 = runif(n, -1, 1),
    exposure = runif(n, 0.6, 3)
  )
  u <- as.vector(as.matrix(d[c("z1", "z2", "z3")]) %*%
                   rep(1 / sqrt(3), 3))
  eta <- 0.5 * d$x1 + sin(2 * u)
  d$y_gaussian <- eta + rnorm(n, sd = 0.5)
  d$y_binomial <- rbinom(n, 1, plogis(eta))
  d$y_poisson <- rpois(n, exp(eta) * d$exposure)
  d
}

fit_phase1_family <- function(fam, data) {
  response <- switch(fam, gaussian = "y_gaussian",
                     binomial = "y_binomial", poisson = "y_poisson")
  args <- list(formula = stats::reformulate("x1", response),
               index = ~ z1 + z2 + z3,
               data = data, family = switch(fam,
                 gaussian = stats::gaussian(), binomial = stats::binomial(),
                 poisson = stats::poisson()),
               M = 5, control = pgplsim_control(maxit = 75, tol = 1e-5))
  if (fam == "poisson") args$offset <- log(data$exposure)
  do.call(pgplsim, args)
}

test_that("all three families retain a coherent public fit/predict/summary/vcov API", {
  d <- make_phase1_compatibility_data()
  for (fam in c("gaussian", "binomial", "poisson")) {
    f <- fit_phase1_family(fam, d)
    expect_s3_class(f, "pgplsim")
    expect_true(isTRUE(f$converged))
    expect_identical(f$family$family, fam)
    expect_true(is.finite(f$lambda) && f$lambda > 0)
    expect_equal(sum(f$alpha_hat^2), 1, tolerance = 1e-6)

    p <- as.numeric(predict(f, type = "response"))
    eta <- as.numeric(predict(f, type = "link"))
    expect_length(p, nrow(d))
    expect_true(all(is.finite(p)))
    if (fam == "gaussian") {
      expect_equal(p, eta, tolerance = 1e-6)
      expect_equal(p, as.numeric(f$fitted_mean), tolerance = 1e-6)
    } else if (fam == "binomial") {
      expect_equal(p, stats::plogis(eta), tolerance = 1e-6)
      expect_true(all(p > 0 & p < 1))
    } else {
      expect_equal(p, exp(eta), tolerance = 1e-6)
      expect_true(all(p > 0))
    }

    ss <- summary(f)
    expect_s3_class(ss, "summary.pgplsim")
    expect_identical(ss$family, fam)
    for (type in c("model", "sandwich", "leverage", "finite")) {
      V <- stats::vcov(f, type = type)
      expect_true(is.matrix(V))
      expect_equal(nrow(V), length(f$alpha_hat) + length(f$beta_hat))
      expect_true(all(is.finite(V)))
      expect_equal(V, t(V), tolerance = 1e-7)
      expect_true(all(diag(V) >= 0))
    }
  }
})

test_that("Gaussian family restrictions and Poisson offset conventions remain explicit", {
  d <- make_phase1_compatibility_data()
  fg <- fit_phase1_family("gaussian", d)
  expect_error(predict(fg, offset = rep(0, nrow(d))), "Poisson")
  expect_error(pgplsim(y_gaussian ~ x1, index = ~ z1 + z2 + z3,
                       data = d, family = gaussian(),
                       offset = rep(0, nrow(d))), "Poisson")
  expect_error(bootstrap_pgplsim(fg, B = 2, method = "stratified"),
               "method = 'case'")
  expect_error(pgplsim(y_gaussian ~ x1, index = ~ z1 + z2 + z3,
                       data = d, family = stats::Gamma()),
               "Family not supported")
  fp <- fit_phase1_family("poisson", d)
  expect_equal(as.numeric(predict(fp, type = "response")),
               as.numeric(fp$fitted_mean), tolerance = 1e-6)
})
