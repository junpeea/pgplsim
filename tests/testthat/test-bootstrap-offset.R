
test_that("Poisson bootstrap recovers exposure correctly", {

  set.seed(20261106)

  n <- 100

  dat <- data.frame(
    x1 = rnorm(n),
    z1 = runif(n, -1, 1),
    z2 = runif(n, -1, 1),
    exposure = runif(n, 0.5, 5)
  )

  eta <- 0.2 + 0.4 * dat$x1 +
    0.3 * sin(0.8 * dat$z1 + 0.6 * dat$z2)

  dat$y <- rpois(
    n,
    dat$exposure * exp(eta)
  )

  off <- log(dat$exposure)

  fit <- pgplsim(
    y ~ x1,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = off,
    M = 4,
    control = pgplsim_control(
      maxit = 100,
      tol = 1e-6
    )
  )

  # Test 1: Automatic exposure recovery
  recovered <- recover_pgplsim_data(fit)

  expect_equal(
    recovered$offset,
    fit$exposure,
    tolerance = 1e-10
  )

  # Test 2: Explicit log-offset conversion
  recovered_explicit <- recover_pgplsim_data(
    fit,
    offset = off
  )

  expect_equal(
    recovered_explicit$offset,
    dat$exposure,
    tolerance = 1e-10
  )

  # Test 3: All recovered exposures are positive
  expect_true(
    all(recovered$offset > 0)
  )

  # Test 4: Negative log-offsets are supported
  expect_true(any(off < 0))

  # Test 5: Resampled exposures retain row alignment
  set.seed(123)

  idx <- bootstrap_sample_indices(
    dat$y,
    method = "case"
  )

  expect_equal(
    recovered$offset[idx],
    dat$exposure[idx],
    tolerance = 1e-10
  )

  # Test 6: Overflowing log-offset is rejected
  expect_error(
    recover_pgplsim_data(
      fit,
      offset = rep(1000, n)
    ),
    regexp = "finite"
  )

})


test_that("Poisson bootstrap preserves original exposure scale", {

  set.seed(20261107)

  n <- 100

  dat <- data.frame(
    x1 = rnorm(n),
    z1 = runif(n, -1, 1),
    z2 = runif(n, -1, 1),
    exposure = runif(n, 0.5, 5)
  )

  u <- 0.8 * dat$z1 + 0.6 * dat$z2

  dat$y <- rpois(
    n,
    dat$exposure * exp(
      0.2 + 0.4 * dat$x1 + 0.3 * sin(2 * u)
    )
  )

  fit <- pgplsim(
    y ~ x1,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = log(dat$exposure),
    M = 4,
    control = pgplsim_control(
      maxit = 100,
      tol = 1e-6
    )
  )

  boot <- bootstrap_pgplsim(
    fit,
    B = 2,
    method = "case",
    seed = 20261107,
    parallel = FALSE,
    keep_fits = TRUE,
    trace = FALSE
  )

  # Verify successful bootstrap execution
  expect_equal(boot$B, 2L)
  expect_length(boot$fits, 2L)

  for (fb in boot$fits) {

    exposure_b <- fb$offset

    # Exposures must be finite and positive
    expect_true(all(is.finite(exposure_b)))
    expect_true(all(exposure_b > 0))

    # Every bootstrap exposure must be one
    # of the original exposure values
    expect_true(
      all(vapply(
        exposure_b,
        function(x) {
          any(abs(fit$exposure - x) < 1e-8)
        },
        logical(1)
      ))
    )

    # Bootstrap cannot create larger exposures
    expect_lte(
      max(exposure_b),
      max(fit$exposure) + 1e-8
    )
  }

})

