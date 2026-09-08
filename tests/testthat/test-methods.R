test_that("print and summary methods run without errors", {

  dat <- make_binomial_test_data()

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  expect_output(
    print(fit),
    "Generalized Partially Linear Single-Index Model"
  )

  fit_summary <- summary(fit)

  expect_s3_class(
    fit_summary,
    "summary.pgplsim"
  )

  expect_output(
    print(fit_summary),
    "Linear component"
  )

  expect_output(
    print(fit_summary),
    "Index component"
  )
})

test_that("plot method runs without error", {

  dat <- make_binomial_test_data()

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  temporary_pdf <- tempfile(
    fileext = ".pdf"
  )

  grDevices::pdf(temporary_pdf)

  expect_no_error(
    plot(fit)
  )

  grDevices::dev.off()

  unlink(temporary_pdf)
})


test_that("summary object has a stable public structure", {

  dat <- make_binomial_test_data(n = 100, seed = 1601)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.pgplsim")

  expect_named(
    s,
    c(
      "call",
      "family",
      "linear",
      "index",
      "lambda",
      "logLik",
      "iterations",
      "converged"
    )
  )

  expect_equal(nrow(s$linear), ncol(dat$X))
  expect_equal(nrow(s$index), ncol(dat$Z))
  expect_true(all(c("Estimate", "SE", "z", "p") %in% names(s$linear)))
  expect_true(all(c("Alpha", "SE") %in% names(s$index)))
  expect_type(s$converged, "logical")
})


test_that("print methods return their inputs invisibly", {

  dat <- make_binomial_test_data(n = 100, seed = 1602)

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  s <- summary(fit)

  expect_invisible(print(fit))
  expect_invisible(print(s))
})
