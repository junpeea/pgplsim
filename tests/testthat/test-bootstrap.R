test_that("bootstrap returns expected estimates", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 789
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot <- bootstrap_pgplsim(
    object = fit,
    B = 3,
    method = "stratified",
    seed = 2026,
    max_attempts = 10,
    trace = FALSE
  )

  expect_s3_class(
    boot,
    "bootstrap.pgplsim"
  )

  expect_equal(
    nrow(boot$estimates_alpha),
    3
  )

  expect_equal(
    ncol(boot$estimates_alpha),
    ncol(dat$Z)
  )

  expect_equal(
    nrow(boot$estimates_beta),
    3
  )

  expect_equal(
    ncol(boot$estimates_beta),
    ncol(dat$X)
  )
})


test_that("bootstrap is reproducible with a fixed seed", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1501
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot1 <- bootstrap_pgplsim(
    object = fit,
    B = 3,
    method = "stratified",
    seed = 2026,
    max_attempts = 12,
    trace = FALSE
  )

  boot2 <- bootstrap_pgplsim(
    object = fit,
    B = 3,
    method = "stratified",
    seed = 2026,
    max_attempts = 12,
    trace = FALSE
  )

  expect_equal(boot1$estimates_alpha, boot2$estimates_alpha)
  expect_equal(boot1$estimates_beta, boot2$estimates_beta)
  expect_equal(boot1$vcov_parameters, boot2$vcov_parameters)
})


test_that("parallel bootstrap returns the expected structure", {

  skip_on_cran()

  cores <- parallel::detectCores(logical = FALSE)
  skip_if(is.na(cores) || cores < 2L, "Two CPU cores are required")

  dat <- make_binomial_test_data(
    n = 80,
    seed = 1502
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot <- bootstrap_pgplsim(
    object = fit,
    B = 2,
    method = "stratified",
    seed = 2026,
    max_attempts = 8,
    parallel = TRUE,
    ncores = 2,
    trace = FALSE
  )

  expect_s3_class(boot, "bootstrap.pgplsim")
  expect_equal(nrow(boot$estimates_alpha), 2)
  expect_equal(nrow(boot$estimates_beta), 2)
})


test_that("bootstrap print method returns its object invisibly", {

  dat <- make_binomial_test_data(
    n = 80,
    seed = 1503
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 25
  )

  boot <- bootstrap_pgplsim(
    object = fit,
    B = 2,
    method = "stratified",
    seed = 2026,
    max_attempts = 8,
    trace = FALSE
  )

  expect_invisible(print(boot))
})
