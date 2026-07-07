############################################################
# Utility functions for LuGPLSIM
############################################################


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
