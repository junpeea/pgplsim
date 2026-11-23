############################################################
# Summary method for pgplsim
############################################################

#' Summarize a pgplsim object
#'
#' @param object An object of class `"pgplsim"`.
#' @param ... Additional arguments, currently ignored.
#'
#' @return An object of class `"summary.pgplsim"`. For Gaussian fits,
#'   includes residual dispersion, residual degrees of freedom, and
#'   asymptotic-normal inference for linear coefficients in the manuscript's
#'   convention.
#'
#' @importFrom stats pnorm
#' @export
summary.pgplsim <- function(object, ...) {

  beta <- object$beta_hat
  se_beta <- object$se_beta

  if (is.null(names(beta))) {
    names(beta) <- paste0("X", seq_along(beta))
  }

  # Lu (2026) establishes asymptotic normality for parametric effects;
  # preserve the package's normal-reference statistics for all families.
  z_beta <- beta / se_beta
  p_beta <- 2 * stats::pnorm(abs(z_beta), lower.tail = FALSE)

  linear_table <- data.frame(
    Estimate = beta,
    SE = se_beta,
    z = z_beta,
    p = p_beta,
    row.names = names(beta),
    check.names = FALSE
  )

  alpha <- object$alpha_hat
  se_alpha <- object$se_alpha

  if (is.null(names(alpha))) {
    names(alpha) <- paste0("Z", seq_along(alpha))
  }

  index_table <- data.frame(
    Alpha = alpha,
    SE = se_alpha,
    row.names = names(alpha),
    check.names = FALSE
  )

  out <- list(
    call = object$call,
    family = object$family$family,
    linear = linear_table,
    index = index_table,
    lambda = object$lambda,
    logLik = object$logLik,
    iterations = object$iterations,
    converged = object$converged
  )

  if (identical(object$family$family, "gaussian")) {
    out$dispersion <- object$dispersion
    out$residual_df <- object$residual_df
    out$effective_df <- object$effective_df
    out$rss <- object$rss
    out$logLik_mle <- object$logLik_mle
  }

  class(out) <- "summary.pgplsim"

  out
}


#' Print summary of a pgplsim object
#'
#' @param x An object of class `"summary.pgplsim"`.
#' @param digits Number of digits to print.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.summary.pgplsim <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...
) {

  cat("\nGeneralized Partially Linear Single-Index Model\n")
  cat("================================================\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nFamily: ", tools::toTitleCase(x$family), "\n", sep = "")
  cat("Log-likelihood: ", signif(x$logLik, digits), "\n", sep = "")
  cat("Smoothing parameter lambda: ", signif(x$lambda, digits), "\n", sep = "")
  cat("Iterations: ", x$iterations, "\n", sep = "")
  cat("Converged: ", x$converged, "\n", sep = "")
  if (identical(x$family, "gaussian")) {
    cat("Residual dispersion: ", signif(x$dispersion, digits), "\n", sep = "")
    cat("Residual df: ", signif(x$residual_df, digits), "\n", sep = "")
    cat("Effective model df: ", signif(x$effective_df, digits), "\n", sep = "")
  }

  cat("\nLinear component:\n")
  print(round(x$linear, digits))

  cat("\nIndex component:\n")
  print(round(x$index, digits))

  invisible(x)
}
