############################################################
# Utility functions for LuGPLSIM
############################################################

alpha_from_theta <- function(theta) {
  denom <- sqrt(1 + sum(theta^2))
  c(1, theta) / denom
}

theta_from_alpha <- function(alpha) {
  stopifnot(length(alpha) >= 2)
  alpha[-1] / alpha[1]
}

jac_alpha_wrt_theta <- function(theta) {
  k <- length(theta)
  s <- sqrt(1 + sum(theta^2))

  J <- matrix(0, nrow = k + 1, ncol = k)
  J[1, ] <- -theta / s^3
  J[2:(k + 1), ] <- diag(1 / s, k) - (theta %o% theta) / s^3

  J
}

var_alpha_from_var_theta <- function(theta, V_theta) {
  J <- jac_alpha_wrt_theta(theta)
  J %*% V_theta %*% t(J)
}

# Normalize single-index coefficient
normalize_alpha <- function(alpha) {

  alpha <- as.numeric(alpha)

  alpha <- alpha / sqrt(sum(alpha^2))

  # sign constraint for identifiability
  if (alpha[1] < 0) {
    alpha <- -alpha
  }

  alpha
}



# Check convergence
check_convergence <- function(old, new, tol = 1e-6) {

  diff <- max(abs(old - new))

  diff < tol
}



# Safe inverse
safe_solve <- function(A, eps = 1e-8) {

  solve(
    A + diag(eps, nrow(A))
  )
}



# Create second difference matrix
difference_matrix <- function(K, order = 2) {

  diff(diag(K), differences = order)
}
