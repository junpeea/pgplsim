test_that("predict returns correct-length binomial predictions", {

  dat <- make_binomial_test_data()

  fit <- pgplsim:::pgplsim_fit(
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

  fit <- pgplsim:::pgplsim_fit(
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


test_that("predict can use stored training matrices", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1301
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  predicted <- predict(
    fit,
    type = "response"
  )

  expect_length(
    predicted,
    length(dat$y)
  )

  expect_equal(
    as.numeric(predicted),
    as.numeric(fit$fitted_mean),
    tolerance = 1e-6
  )
})


test_that("predict validates new-data dimensions", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1302
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  expect_error(
    predict(
      fit,
      newX = dat$X[, 1, drop = FALSE],
      newZ = dat$Z
    ),
    "linear prediction matrix"
  )

  expect_error(
    predict(
      fit,
      newX = dat$X,
      newZ = dat$Z[, 1:2, drop = FALSE]
    ),
    "index prediction matrix"
  )

  expect_error(
    predict(
      fit,
      newX = dat$X[-1, , drop = FALSE],
      newZ = dat$Z
    ),
    "same number of rows"
  )
})


test_that("predict rejects missing and non-finite new data", {

  dat <- make_binomial_test_data(
    n = 100,
    seed = 1303
  )

  fit <- pgplsim:::pgplsim_fit(
    y = dat$y,
    X = dat$X,
    Z = dat$Z,
    family = binomial(),
    M = 4,
    maxit = 30
  )

  bad_X <- dat$X
  bad_X[1, 1] <- NA_real_

  expect_error(
    predict(
      fit,
      newX = bad_X,
      newZ = dat$Z
    ),
    "finite, non-missing"
  )

  bad_Z <- dat$Z
  bad_Z[1, 1] <- Inf

  expect_error(
    predict(
      fit,
      newX = dat$X,
      newZ = bad_Z
    ),
    "finite, non-missing"
  )
})


test_that("Poisson prediction validates offset length and values", {

  set.seed(1304)

  n <- 120

  X <- cbind(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  Z <- cbind(
    Z1 = rnorm(n),
    Z2 = rnorm(n),
    Z3 = rnorm(n)
  )

  exposure <- runif(
    n,
    min = 100,
    max = 1000
  )

  eta <- -5 +
    0.4 * X[, "X1"] -
    0.2 * X[, "X2"] +
    0.3 * sin(
      Z[, "Z1"] +
        Z[, "Z2"]
    )

  y <- rpois(
    n,
    lambda = exposure * exp(eta)
  )

  fit <- pgplsim:::pgplsim_fit(
    y = y,
    X = X,
    Z = Z,
    family = poisson(),

    # Internal engine uses positive exposure scale.
    offset = exposure,

    M = 4,
    maxit = 30
  )

  ############################################################
  # Incorrect offset length
  ############################################################

  expect_error(
    predict(
      fit,
      newX = X,
      newZ = Z,
      offset = rep(0, n - 1)
    ),
    "`offset` must have"
  )

  ############################################################
  # Missing offset value
  ############################################################

  bad_offset <- rep(0, n)
  bad_offset[1] <- NA_real_

  expect_error(
    predict(
      fit,
      newX = X,
      newZ = Z,
      offset = bad_offset
    ),
    "finite, non-missing"
  )

  ############################################################
  # Infinite offset value
  ############################################################

  bad_offset <- rep(0, n)
  bad_offset[1] <- Inf

  expect_error(
    predict(
      fit,
      newX = X,
      newZ = Z,
      offset = bad_offset
    ),
    "finite, non-missing"
  )

  ############################################################
  # Negative link-scale offsets are valid
  ############################################################

  expect_no_error(
    predict(
      fit,
      newX = X,
      newZ = Z,
      offset = rep(-2, n)
    )
  )

  ############################################################
  # Zero link-scale offset is valid
  ############################################################

  expect_no_error(
    predict(
      fit,
      newX = X,
      newZ = Z,
      offset = rep(0, n)
    )
  )
})

make_formula_binomial_prediction_data <- function(
    n = 200,
    seed = 2026
) {

  set.seed(seed)

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    group = factor(
      sample(
        c("A", "B", "C"),
        n,
        replace = TRUE
      ),
      levels = c("A", "B", "C")
    ),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  group_effect <- c(
    A = 0,
    B = 0.3,
    C = -0.2
  )

  eta <- -0.2 +
    0.5 * dat$x1 -
    0.3 * dat$x2 +
    group_effect[dat$group] +
    0.4 * sin(dat$z1 + dat$z2)

  dat$y <- rbinom(
    n,
    size = 1,
    prob = plogis(eta)
  )

  dat
}

test_that("formula-based binomial models predict from newdata", {

  dat <- make_formula_binomial_prediction_data()

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4,
    control = pgplsim_control(
      maxit = 30,
      tol = 1e-5
    )
  )

  new_dat <- dat[1:25, , drop = FALSE]

  predicted_response <- predict(
    fit,
    newdata = new_dat,
    type = "response"
  )

  predicted_link <- predict(
    fit,
    newdata = new_dat,
    type = "link"
  )

  predicted_index <- predict(
    fit,
    newdata = new_dat,
    type = "index"
  )

  predicted_smooth <- predict(
    fit,
    newdata = new_dat,
    type = "smooth"
  )

  predicted_linear <- predict(
    fit,
    newdata = new_dat,
    type = "linear"
  )

  expect_length(
    predicted_response,
    nrow(new_dat)
  )

  expect_length(
    predicted_link,
    nrow(new_dat)
  )

  expect_length(
    predicted_index,
    nrow(new_dat)
  )

  expect_length(
    predicted_smooth,
    nrow(new_dat)
  )

  expect_length(
    predicted_linear,
    nrow(new_dat)
  )

  expect_true(
    all(predicted_response >= 0)
  )

  expect_true(
    all(predicted_response <= 1)
  )

  expect_equal(
    predicted_link,
    predicted_linear + predicted_smooth,
    tolerance = 1e-8
  )
})

test_that("formula newdata prediction agrees with manual matrices", {

  dat <- make_formula_binomial_prediction_data(
    n = 180,
    seed = 2027
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  new_dat <- dat[31:60, , drop = FALSE]

  ############################################################
  # Public formula-based prediction
  ############################################################

  pred_newdata <- predict(
    fit,
    newdata = new_dat,
    type = "response"
  )

  ############################################################
  # Equivalent manually constructed matrices
  ############################################################

  newX <- model.matrix(
    ~ x1 + x2,
    data = new_dat
  )

  newZ <- model.matrix(
    ~ z1 + z2 - 1,
    data = new_dat
  )

  pred_matrix <- predict(
    fit,
    newX = newX,
    newZ = newZ,
    type = "response"
  )

  expect_equal(
    pred_newdata,
    pred_matrix,
    tolerance = 1e-8
  )
})

test_that("formula prediction handles factor predictors", {

  dat <- make_formula_binomial_prediction_data(
    n = 240,
    seed = 2028
  )

  fit <- pgplsim(
    y ~ x1 + x2 + group,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  ############################################################
  # Include observations from all fitted factor levels
  ############################################################

  idx_A <- which(dat$group == "A")[1]
  idx_B <- which(dat$group == "B")[1]
  idx_C <- which(dat$group == "C")[1]

  new_dat <- dat[
    c(idx_A, idx_B, idx_C),
    ,
    drop = FALSE
  ]

  pred_newdata <- predict(
    fit,
    newdata = new_dat,
    type = "response"
  )

  ############################################################
  # Reconstruct expected linear design matrix
  ############################################################

  newX <- model.matrix(
    ~ x1 + x2 + group,
    data = new_dat,
    contrasts.arg = fit$contrasts
  )

  newZ <- model.matrix(
    ~ z1 + z2 - 1,
    data = new_dat
  )

  pred_matrix <- predict(
    fit,
    newX = newX,
    newZ = newZ,
    type = "response"
  )

  expect_length(
    pred_newdata,
    3
  )

  expect_true(
    all(is.finite(pred_newdata))
  )

  expect_equal(
    pred_newdata,
    pred_matrix,
    tolerance = 1e-8
  )
})

test_that("formula-based Poisson prediction uses link-scale offset", {

  set.seed(2029)

  n <- 220

  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n),
    exposure = runif(
      n,
      min = 100,
      max = 1000
    )
  )

  eta <- -5 +
    0.4 * dat$x1 -
    0.2 * dat$x2 +
    0.4 * sin(dat$z1 + dat$z2)

  dat$y <- rpois(
    n,
    lambda = dat$exposure * exp(eta)
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = poisson(),
    offset = log(dat$exposure),
    M = 4,
    control = pgplsim_control(
      maxit = 30
    )
  )

  new_dat <- dat[1:30, , drop = FALSE]

  ############################################################
  # Public prediction using link-scale offset
  ############################################################

  pred_response <- predict(
    fit,
    newdata = new_dat,
    offset = log(new_dat$exposure),
    type = "response"
  )

  pred_link <- predict(
    fit,
    newdata = new_dat,
    offset = log(new_dat$exposure),
    type = "link"
  )

  ############################################################
  # Response and link predictions must agree
  ############################################################

  expect_equal(
    pred_response,
    exp(pred_link),
    tolerance = 1e-8
  )

  expect_true(
    all(pred_response > 0)
  )

  ############################################################
  # Removing exposure corresponds to offset = 0
  ############################################################

  pred_unit_exposure <- predict(
    fit,
    newdata = new_dat,
    offset = 0,
    type = "response"
  )

  expect_equal(
    pred_response,
    pred_unit_exposure * new_dat$exposure,
    tolerance = 1e-6
  )
})

test_that("predict rejects newdata together with newX or newZ", {

  dat <- make_formula_binomial_prediction_data(
    n = 120,
    seed = 2030
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  new_dat <- dat[1:10, , drop = FALSE]

  newX <- model.matrix(
    ~ x1 + x2,
    data = new_dat
  )

  newZ <- model.matrix(
    ~ z1 + z2 - 1,
    data = new_dat
  )

  expect_error(
    predict(
      fit,
      newdata = new_dat,
      newX = newX
    ),
    "either `newdata` or `newX`/`newZ`"
  )

  expect_error(
    predict(
      fit,
      newdata = new_dat,
      newZ = newZ
    ),
    "either `newdata` or `newX`/`newZ`"
  )
})

test_that("formula prediction rejects missing variables in newdata", {

  dat <- make_formula_binomial_prediction_data(
    n = 120,
    seed = 2031
  )

  fit <- pgplsim(
    y ~ x1 + x2,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  ############################################################
  # Missing linear predictor
  ############################################################

  new_dat <- dat[1:10, , drop = FALSE]
  new_dat$x2 <- NULL

  expect_error(
    predict(
      fit,
      newdata = new_dat
    ),
    "linear prediction data"
  )

  ############################################################
  # Missing index predictor
  ############################################################

  new_dat <- dat[1:10, , drop = FALSE]
  new_dat$z2 <- NULL

  expect_error(
    predict(
      fit,
      newdata = new_dat
    ),
    "index prediction data"
  )
})

test_that("formula prediction rejects unseen factor levels", {

  dat <- make_formula_binomial_prediction_data(
    n = 180,
    seed = 2032
  )

  fit <- pgplsim(
    y ~ x1 + group,
    index = ~ z1 + z2,
    data = dat,
    family = binomial(),
    M = 4
  )

  new_dat <- dat[1:10, , drop = FALSE]

  new_dat$group <- as.character(
    new_dat$group
  )

  new_dat$group[1] <- "D"

  new_dat$group <- factor(
    new_dat$group
  )

  expect_error(
    predict(
      fit,
      newdata = new_dat
    ),
    "new level"
  )
})
