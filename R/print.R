############################################################
# Print method for LuGPLSIM
############################################################

#' Print a LuGPLSIM object
#'
#' @param x An object of class `"LuGPLSIM"`.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.LuGPLSIM <- function(x, ...) {

  fam <- if (!is.null(x$family$family)) x$family$family else "unknown"

  cat("\nGeneralized Partially Linear Single-Index Model\n")
  cat("------------------------------------------------\n")
  cat("Family:          ", tools::toTitleCase(fam), "\n", sep = "")
  cat("Index variables: ", length(x$alpha_hat), "\n", sep = "")
  cat("Linear variables:", length(x$beta_hat), "\n", sep = " ")

  cat("\nSmoothing parameter:\n")
  cat("lambda = ", signif(x$lambda, 4), "\n", sep = "")

  cat("\nModel fit:\n")
  cat("Iterations: ", x$iterations, "\n", sep = "")
  cat("Converged:  ", x$converged, "\n", sep = "")

  invisible(x)
}
