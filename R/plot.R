############################################################
# Plot method for LuGPLSIM
############################################################

#' Plot a fitted LuGPLSIM model
#'
#' Produces a plot of the estimated smooth single-index component
#' against the fitted index values.
#'
#' @param x An object of class `"LuGPLSIM"`.
#' @param type Type of plot. Currently `"smooth"` is supported.
#' @param rug Logical; if `TRUE`, add a rug plot of the fitted index values.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional plot title.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns the plotted data.
#'
#' @export
plot.LuGPLSIM <- function(
    x,
    type = c("smooth"),
    rug = TRUE,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ...
) {

  if (!inherits(x, "LuGPLSIM")) {
    stop(
      "x must inherit from class \"LuGPLSIM\".",
      call. = FALSE
    )
  }

  type <- match.arg(type)

  if (is.null(x$u_hat)) {
    stop(
      "The fitted object does not contain u_hat.",
      call. = FALSE
    )
  }

  u_hat <- as.numeric(x$u_hat)

  if (anyNA(u_hat) || any(!is.finite(u_hat))) {
    stop(
      "u_hat contains non-finite or missing values.",
      call. = FALSE
    )
  }

  ############################################################
  # Recover the estimated smooth component
  ############################################################

  smooth_hat <- NULL

  # Preferred: directly stored smooth component
  if (
    !is.null(x$smooth_component) &&
    length(x$smooth_component) == length(u_hat)
  ) {
    smooth_hat <- as.numeric(x$smooth_component)
  }

  # Second choice: transformed fitted basis matrix
  if (
    is.null(smooth_hat) &&
    !is.null(x$basis_matrix)
  ) {

    basis_matrix <- as.matrix(x$basis_matrix)

    if (
      nrow(basis_matrix) == length(u_hat) &&
      ncol(basis_matrix) == length(x$gamma_hat)
    ) {
      smooth_hat <- as.numeric(
        basis_matrix %*% x$gamma_hat
      )
    }
  }

  # Backward-compatible fitted basis name
  if (
    is.null(smooth_hat) &&
    !is.null(x$BB_hat)
  ) {

    BB_hat <- as.matrix(x$BB_hat)

    if (
      nrow(BB_hat) == length(u_hat) &&
      ncol(BB_hat) == length(x$gamma_hat)
    ) {
      smooth_hat <- as.numeric(
        BB_hat %*% x$gamma_hat
      )
    }
  }

  # Final fallback: reconstruct from fitted linear predictor
  if (
    is.null(smooth_hat) &&
    !is.null(x$linear_predictors) &&
    !is.null(x$X) &&
    !is.null(x$beta_hat)
  ) {

    X <- as.matrix(x$X)

    if (
      nrow(X) == length(u_hat) &&
      ncol(X) == length(x$beta_hat) &&
      length(x$linear_predictors) == length(u_hat)
    ) {
      smooth_hat <- as.numeric(
        x$linear_predictors -
          X %*% x$beta_hat
      )
    }
  }

  if (is.null(smooth_hat)) {
    stop(
      paste0(
        "The estimated smooth component could not be recovered. ",
        "The fitted object should contain smooth_component, ",
        "basis_matrix, BB_hat, or compatible linear_predictors."
      ),
      call. = FALSE
    )
  }

  if (
    anyNA(smooth_hat) ||
    any(!is.finite(smooth_hat))
  ) {
    stop(
      "The estimated smooth component contains non-finite values.",
      call. = FALSE
    )
  }

  ############################################################
  # Order by fitted index for plotting
  ############################################################

  ordering <- order(u_hat)

  plot_data <- data.frame(
    index = u_hat[ordering],
    smooth = smooth_hat[ordering]
  )

  if (is.null(xlab)) {
    xlab <- expression(
      hat(u) == Z^T * hat(alpha)
    )
  }

  if (is.null(ylab)) {
    ylab <- expression(
      hat(phi)(hat(u))
    )
  }

  if (is.null(main)) {
    main <- "Estimated single-index smooth"
  }

  graphics::plot(
    plot_data$index,
    plot_data$smooth,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  if (isTRUE(rug)) {
    graphics::rug(u_hat)
  }

  invisible(plot_data)
}
