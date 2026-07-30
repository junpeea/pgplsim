make_binomial_test_data <- function(
    n = 150,
    seed = 123
) {
  set.seed(seed)

  Z1 <- runif(n, 0, 2)
  Z2 <- runif(n, 0, 2)
  Z3 <- runif(n, 0, 2)

  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)

  Z <- cbind(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3
  )

  X <- cbind(
    X1 = X1,
    X2 = X2
  )

  alpha_true <- c(1, 1, 1) / sqrt(3)
  beta_true <- c(1, -0.5)

  u <- as.numeric(Z %*% alpha_true)

  phi <- sin(2 * u)
  phi <- phi - mean(phi)

  linear_predictor <- as.numeric(
    X %*% beta_true + phi
  )

  probability <- plogis(linear_predictor)

  y <- rbinom(
    n = n,
    size = 1,
    prob = probability
  )

  list(
    y = y,
    X = X,
    Z = Z,
    alpha_true = alpha_true,
    beta_true = beta_true
  )
}


make_poisson_test_data <- function(
    n = 150,
    seed = 456
) {
  set.seed(seed)

  Z1 <- runif(n, 0, 2)
  Z2 <- runif(n, 0, 2)
  Z3 <- runif(n, 0, 2)

  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)

  exposure <- runif(n, 1, 5)

  Z <- cbind(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3
  )

  X <- cbind(
    X1 = X1,
    X2 = X2
  )

  alpha_true <- c(1, 1, 1) / sqrt(3)
  beta_true <- c(0.4, -0.3)

  u <- as.numeric(Z %*% alpha_true)

  phi <- 0.5 * sin(2 * u)
  phi <- phi - mean(phi)

  linear_predictor <- as.numeric(
    X %*% beta_true + phi
  )

  mu <- exposure * exp(linear_predictor)

  y <- rpois(
    n = n,
    lambda = mu
  )

  list(
    y = y,
    X = X,
    Z = Z,
    exposure = exposure,
    alpha_true = alpha_true,
    beta_true = beta_true
  )
}
