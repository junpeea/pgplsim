#' Predict from a fitted LuGPLSIM model
#'
#' @param object A fitted object of class `"LuGPLSIM"`.
#' @param newX Optional matrix of new linear covariates. If omitted,
#'   the original fitting matrix is used.
#' @param newZ Optional matrix of new single-index covariates. If omitted,
#'   the original fitting matrix is used.
#' @param type Type of prediction. One of `"response"`, `"link"`,
#'   `"index"`, `"smooth"`, or `"linear"`.
#' @param offset Optional exposure vector for Poisson prediction.
#'   This is an exposure on the original scale, not a log offset.
#'   If omitted for in-sample prediction, the fitted exposure is used.
#'   If omitted for new Poisson data, exposure is set to one.
#' @param ... Additional arguments, currently unused.
#'
#' @return A numeric vector of predictions.
#'
#' @importFrom splines splineDesign
#' @export
predict.LuGPLSIM <- function(
    object,
    newX = NULL,
    newZ = NULL,
    type = c(
      "response",
      "link",
      "index",
      "smooth",
      "linear"
    ),
    offset = NULL,
    ...
) {

  type <- match.arg(type)

  ############################################################
  # 1. Validate fitted object
  ############################################################

  if (!inherits(object, "LuGPLSIM")) {
    stop(
      "object must inherit from class \"LuGPLSIM\".",
      call. = FALSE
    )
  }

  required_components <- c(
    "alpha_hat",
    "beta_hat",
    "gamma_hat",
    "family"
  )

  missing_components <- required_components[
    vapply(
      required_components,
      function(component_name) {
        is.null(object[[component_name]])
      },
      logical(1)
    )
  ]

  if (length(missing_components) > 0L) {
    stop(
      paste0(
        "The fitted object is missing required component(s): ",
        paste(
          missing_components,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # 2. Obtain prediction matrices
  ############################################################

  if (is.null(newX)) {
    if (is.null(object$X)) {
      stop(
        paste0(
          "newX was not supplied and the fitted object ",
          "does not contain the original X matrix."
        ),
        call. = FALSE
      )
    }

    newX <- object$X
  }

  if (is.null(newZ)) {
    if (is.null(object$Z)) {
      stop(
        paste0(
          "newZ was not supplied and the fitted object ",
          "does not contain the original Z matrix."
        ),
        call. = FALSE
      )
    }

    newZ <- object$Z
  }

  newX <- as.matrix(newX)
  newZ <- as.matrix(newZ)

  if (nrow(newX) != nrow(newZ)) {
    stop(
      "newX and newZ must have the same number of rows.",
      call. = FALSE
    )
  }

  if (ncol(newX) != length(object$beta_hat)) {
    stop(
      paste0(
        "newX has ",
        ncol(newX),
        " columns, but the fitted model requires ",
        length(object$beta_hat),
        "."
      ),
      call. = FALSE
    )
  }

  if (ncol(newZ) != length(object$alpha_hat)) {
    stop(
      paste0(
        "newZ has ",
        ncol(newZ),
        " columns, but the fitted model requires ",
        length(object$alpha_hat),
        "."
      ),
      call. = FALSE
    )
  }

  if (
    anyNA(newX) ||
    anyNA(newZ) ||
    any(!is.finite(newX)) ||
    any(!is.finite(newZ))
  ) {
    stop(
      paste0(
        "newX and newZ must contain only finite, ",
        "non-missing values."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # 3. Determine whether prediction data equal training data
  ############################################################

  same_X <- !is.null(object$X) &&
    identical(
      dim(newX),
      dim(as.matrix(object$X))
    ) &&
    isTRUE(
      all.equal(
        unname(newX),
        unname(as.matrix(object$X)),
        check.attributes = FALSE
      )
    )

  same_Z <- !is.null(object$Z) &&
    identical(
      dim(newZ),
      dim(as.matrix(object$Z))
    ) &&
    isTRUE(
      all.equal(
        unname(newZ),
        unname(as.matrix(object$Z)),
        check.attributes = FALSE
      )
    )

  in_sample <- same_X && same_Z

  ############################################################
  # 4. Compute single-index values
  ############################################################

  index <- as.numeric(
    newZ %*% object$alpha_hat
  )

  if (identical(type, "index")) {
    return(index)
  }

  ############################################################
  # 5. Reconstruct transformed spline basis
  ############################################################

  BB_new <- NULL

  # For exact in-sample prediction, use the basis actually used
  # in the final fitted model whenever available.
  if (
    in_sample &&
    !is.null(object$basis_matrix)
  ) {
    candidate_basis <- as.matrix(
      object$basis_matrix
    )

    if (
      nrow(candidate_basis) == nrow(newZ) &&
      ncol(candidate_basis) ==
      length(object$gamma_hat)
    ) {
      BB_new <- candidate_basis
    }
  }

  # Backward-compatible in-sample basis name.
  if (
    is.null(BB_new) &&
    in_sample &&
    !is.null(object$BB_hat)
  ) {
    candidate_basis <- as.matrix(
      object$BB_hat
    )

    if (
      nrow(candidate_basis) == nrow(newZ) &&
      ncol(candidate_basis) ==
      length(object$gamma_hat)
    ) {
      BB_new <- candidate_basis
    }
  }

  # For new data, reconstruct raw splineDesign basis and then
  # apply the original Q2 centering transformation.
  if (is.null(BB_new)) {

    basis_spec <- object$basis_spec

    # Allow older objects that stored the relevant components
    # separately but did not yet contain basis_spec.
    if (is.null(basis_spec)) {

      if (
        !is.null(object$knots) &&
        !is.null(object$Q2) &&
        !is.null(object$ord) &&
        !is.null(object$outer_ok)
      ) {
        basis_spec <- list(
          basis_function = "splineDesign",
          knots = object$knots,
          ord = object$ord,
          outer_ok = object$outer_ok,
          Q2 = object$Q2,
          lower_boundary = min(object$knots),
          upper_boundary = max(object$knots)
        )
      } else {
        stop(
          paste0(
            "The fitted object does not contain the spline ",
            "specification required for prediction. Refit the ",
            "model using the current PSPLINE version."
          ),
          call. = FALSE
        )
      }
    }

    required_spec <- c(
      "knots",
      "ord",
      "outer_ok",
      "Q2"
    )

    missing_spec <- required_spec[
      !required_spec %in% names(basis_spec)
    ]

    if (length(missing_spec) > 0L) {
      stop(
        paste0(
          "object$basis_spec is incomplete. Missing: ",
          paste(
            missing_spec,
            collapse = ", "
          ),
          "."
        ),
        call. = FALSE
      )
    }

    knots <- as.numeric(
      basis_spec$knots
    )

    Q2 <- as.matrix(
      basis_spec$Q2
    )

    ord <- as.integer(
      basis_spec$ord
    )

    outer_ok <- isTRUE(
      basis_spec$outer_ok
    )

    lower_boundary <- min(knots)
    upper_boundary <- max(knots)

    if (
      !outer_ok &&
      any(
        index < lower_boundary |
        index > upper_boundary
      )
    ) {
      stop(
        paste0(
          "Some new single-index values lie outside the fitted ",
          "spline boundary [",
          signif(lower_boundary, 6),
          ", ",
          signif(upper_boundary, 6),
          "]. Extrapolation is unavailable because the fitted ",
          "basis used outer.ok = FALSE."
        ),
        call. = FALSE
      )
    }

    B_new <- tryCatch(
      splines::splineDesign(
        knots = knots,
        x = index,
        ord = ord,
        outer.ok = outer_ok
      ),
      error = function(e) {
        stop(
          paste0(
            "The prediction spline basis could not be ",
            "constructed: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    B_new <- as.matrix(B_new)

    if (ncol(B_new) != nrow(Q2)) {
      stop(
        paste0(
          "The raw prediction basis has ",
          ncol(B_new),
          " columns, but Q2 has ",
          nrow(Q2),
          " rows. The stored spline specification is ",
          "internally inconsistent."
        ),
        call. = FALSE
      )
    }

    BB_new <- B_new %*% Q2
  }

  BB_new <- as.matrix(BB_new)

  if (nrow(BB_new) != nrow(newZ)) {
    stop(
      paste0(
        "The transformed prediction basis has ",
        nrow(BB_new),
        " rows, but newZ has ",
        nrow(newZ),
        " rows."
      ),
      call. = FALSE
    )
  }

  if (
    ncol(BB_new) !=
    length(object$gamma_hat)
  ) {
    stop(
      paste0(
        "The transformed prediction basis has ",
        ncol(BB_new),
        " columns, whereas gamma_hat has length ",
        length(object$gamma_hat),
        "."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # 6. Compute model components
  ############################################################

  smooth_component <- as.numeric(
    BB_new %*% object$gamma_hat
  )

  if (identical(type, "smooth")) {
    return(smooth_component)
  }

  linear_component <- as.numeric(
    newX %*% object$beta_hat
  )

  if (identical(type, "linear")) {
    return(linear_component)
  }

  eta <- linear_component +
    smooth_component

  if (identical(type, "link")) {
    return(eta)
  }

  ############################################################
  # 7. Response-scale prediction
  ############################################################

  family_name <- tolower(
    as.character(
      object$family$family
    )[1L]
  )

  if (identical(family_name, "binomial")) {
    return(
      stats::plogis(eta)
    )
  }

  if (identical(family_name, "poisson")) {

    n_new <- nrow(newX)

    if (is.null(offset)) {

      if (
        in_sample &&
        !is.null(object$offset) &&
        length(object$offset) == n_new
      ) {
        exposure <- as.numeric(
          object$offset
        )
      } else {
        exposure <- rep(1, n_new)
      }

    } else {

      exposure <- as.numeric(offset)

      if (length(exposure) == 1L) {
        exposure <- rep(
          exposure,
          n_new
        )
      }

      if (length(exposure) != n_new) {
        stop(
          paste0(
            "For Poisson prediction, offset/exposure must ",
            "have length 1 or the same number of rows as newX."
          ),
          call. = FALSE
        )
      }
    }

    if (
      anyNA(exposure) ||
      any(!is.finite(exposure)) ||
      any(exposure <= 0)
    ) {
      stop(
        paste0(
          "For Poisson prediction, offset/exposure values ",
          "must be finite and strictly positive."
        ),
        call. = FALSE
      )
    }

    eta_bounded <- pmin(
      pmax(eta, -30),
      30
    )

    return(
      exposure * exp(eta_bounded)
    )
  }

  if (identical(family_name, "gaussian")) {
    return(eta)
  }

  if (
    !is.null(object$family$linkinv) &&
    is.function(object$family$linkinv)
  ) {
    return(
      object$family$linkinv(eta)
    )
  }

  stop(
    paste0(
      "Response-scale prediction is not implemented for family: ",
      family_name,
      "."
    ),
    call. = FALSE
  )
}
