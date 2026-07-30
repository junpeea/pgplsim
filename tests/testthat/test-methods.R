test_that("print and summary methods run without errors", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
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
    "summary.LuGPLSIM"
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
