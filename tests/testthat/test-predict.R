test_that("predict returns correct-length binomial predictions", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  predicted_response <- predict(
    fit,
    newX = dat$X,
    newZ = dat$Z,
    type = "response"
  )

  predicted_link <- predict(
    fit,
    newX = dat$X,
    newZ = dat$Z,
    type = "link"
  )

  predicted_index <- predict(
    fit,
    newX = dat$X,
    newZ = dat$Z,
    type = "index"
  )

  expect_length(
    predicted_response,
    length(dat$y)
  )

  expect_length(
    predicted_link,
    length(dat$y)
  )

  expect_length(
    predicted_index,
    length(dat$y)
  )

  expect_true(
    all(predicted_response >= 0)
  )

  expect_true(
    all(predicted_response <= 1)
  )
})

test_that("in-sample predictions match stored fitted means", {

  dat <- make_binomial_test_data()

  fit <- LuGPLSIM(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  predicted <- predict(
    fit,
    newX = dat$X,
    newZ = dat$Z,
    type = "response"
  )

  expect_equal(
    as.numeric(predicted),
    as.numeric(fit$fitted_mean),
    tolerance = 1e-6
  )
})
