#' Predict from a fitted pgplsim model
#'
#' Generate predictions from a fitted generalized partially linear
#' single-index model.
#'
#' For models fitted with the formula interface, new observations should
#' normally be supplied through `newdata`. The linear and single-index design
#' matrices are then reconstructed automatically from the formulas stored in
#' the fitted model.
#'
#' The `newX` and `newZ` arguments provide a low-level matrix interface.
#' For ordinary use, `newdata` is recommended. The matrix arguments should
#' not be supplied together with `newdata`.
#'
#' @param object A fitted object of class `"pgplsim"`.
#' @param newdata Optional data frame containing variables used in the original
#'   model. If `NULL`, prediction uses `newX`/`newZ` when supplied; otherwise it
#'   uses the original design matrices stored in `object`.
#' @param newX Optional matrix of linear-component covariates for low-level
#'   matrix-based prediction. Do not supply together with `newdata`.
#' @param newZ Optional matrix of single-index covariates for low-level
#'   matrix-based prediction. Do not supply together with `newdata`.
#' @param type Character string specifying the prediction scale. One of
#'   `"response"`, `"link"`, `"index"`, `"smooth"`, or `"linear"`.
#' @param offset Optional numeric offset on the linear-predictor scale for
#'   Poisson prediction. A scalar is recycled to all prediction observations.
#'   For an exposure such as population or person-time, supply its logarithm,
#'   for example `offset = log(exposure)`. Offsets are not supported for
#'   binomial models.
#' @param ... Additional arguments passed to the internal matrix prediction
#'   engine.
#'
#' @return A numeric vector of predictions. For `newdata` predictions, rows
#'   containing missing required covariates (or a missing Poisson offset) are
#'   returned as `NA` while preserving the row order of `newdata`.
#'
#' @export
predict.pgplsim <- function(
    object,
    newdata = NULL,
    newX = NULL,
    newZ = NULL,
    type = c("response", "link", "index", "smooth", "linear"),
    offset = NULL,
    ...
) {

  type <- match.arg(type)

  ############################################################
  # 1. Validate fitted object and interface choice
  ############################################################

  if (!inherits(object, "pgplsim")) {
    stop(
      "`object` must inherit from class \"pgplsim\".",
      call. = FALSE
    )
  }

  if (is.null(object$family) || is.null(object$family$family)) {
    stop(
      "The fitted object does not contain a valid family specification.",
      call. = FALSE
    )
  }

  fam <- tolower(as.character(object$family$family)[1L])

  if (!fam %in% c("binomial", "poisson")) {
    stop(
      "Prediction is currently supported only for binomial and Poisson models.",
      call. = FALSE
    )
  }

  if (!is.null(newdata) && (!is.null(newX) || !is.null(newZ))) {
    stop(
      "Supply either `newdata` or `newX`/`newZ`, not both.",
      call. = FALSE
    )
  }

  ############################################################
  # 2. Matrix-based / in-sample prediction
  ############################################################

  if (is.null(newdata)) {

    if (fam == "binomial" && !is.null(offset)) {
      stop(
        "`offset` is not supported for binomial models.",
        call. = FALSE
      )
    }

    return(
      .predict_pgplsim_matrix(
        object = object,
        newX = newX,
        newZ = newZ,
        type = type,
        offset = offset,
        ...
      )
    )
  }

  ############################################################
  # 3. Validate newdata
  ############################################################

  if (!is.data.frame(newdata)) {
    stop(
      "`newdata` must be a data frame.",
      call. = FALSE
    )
  }

  if (nrow(newdata) < 1L) {
    stop(
      "`newdata` must contain at least one observation.",
      call. = FALSE
    )
  }

  if (is.null(object$terms) || is.null(object$index_terms)) {
    stop(
      paste0(
        "`newdata` prediction requires a model fitted with the formula ",
        "interface. For older matrix-based fits, supply `newX` and `newZ` ",
        "instead."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # 4. Reconstruct linear-component model frame
  ############################################################

  linear_terms <- stats::delete.response(object$terms)

  mf_linear <- tryCatch(
    stats::model.frame(
      formula = linear_terms,
      data = newdata,
      na.action = stats::na.pass,
      xlev = object$xlevels,
      drop.unused.levels = FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct the linear prediction data: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  ############################################################
  # 5. Reconstruct single-index model frame
  ############################################################

  mf_index <- tryCatch(
    stats::model.frame(
      formula = object$index_terms,
      data = newdata,
      na.action = stats::na.pass,
      xlev = object$index_xlevels,
      drop.unused.levels = FALSE
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct the index prediction data: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  ############################################################
  # 6. Common complete-case prediction sample
  ############################################################

  complete <- stats::complete.cases(mf_linear) &
    stats::complete.cases(mf_index)

  ############################################################
  # 7. Validate prediction offset
  ############################################################

  if (fam == "binomial") {

    if (!is.null(offset)) {
      stop(
        "`offset` is not supported for binomial models.",
        call. = FALSE
      )
    }

  } else if (!is.null(offset)) {

    if (!is.numeric(offset)) {
      stop(
        "`offset` must be numeric.",
        call. = FALSE
      )
    }

    if (!length(offset) %in% c(1L, nrow(newdata))) {
      stop(
        paste0(
          "For Poisson prediction, `offset` must have length 1 or ",
          "the same number of observations as `newdata`."
        ),
        call. = FALSE
      )
    }

    if (length(offset) == 1L) {
      offset <- rep(offset, nrow(newdata))
    }

    if (any(is.infinite(offset))) {
      stop(
        "`offset` must not contain infinite values.",
        call. = FALSE
      )
    }

    complete <- complete & !is.na(offset)
  }

  ############################################################
  # 8. Return NA if all newdata rows are incomplete
  ############################################################

  if (!any(complete)) {
    out <- rep(NA_real_, nrow(newdata))
    return(out)
  }

  mf_linear_complete <- mf_linear[complete, , drop = FALSE]
  mf_index_complete <- mf_index[complete, , drop = FALSE]

  ############################################################
  # 9. Construct linear prediction matrix X
  ############################################################

  newX <- tryCatch(
    stats::model.matrix(
      object = linear_terms,
      data = mf_linear_complete,
      contrasts.arg = object$contrasts
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct the linear prediction matrix: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  ############################################################
  # 10. Construct single-index prediction matrix Z
  ############################################################

  newZ <- tryCatch(
    stats::model.matrix(
      object = object$index_terms,
      data = mf_index_complete,
      contrasts.arg = object$index_contrasts
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct the index prediction matrix: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  ############################################################
  # 11. Defensive design-column checks
  ############################################################

  expected_X <- colnames(object$X)
  expected_Z <- colnames(object$Z)

  if (!is.null(expected_X)) {

    missing_X <- setdiff(expected_X, colnames(newX))
    extra_X <- setdiff(colnames(newX), expected_X)

    if (length(missing_X) > 0L || length(extra_X) > 0L) {
      stop(
        "The linear-component design matrix generated from `newdata` does not match the fitted model.",
        call. = FALSE
      )
    }

    newX <- newX[, expected_X, drop = FALSE]
  }

  if (!is.null(expected_Z)) {

    missing_Z <- setdiff(expected_Z, colnames(newZ))
    extra_Z <- setdiff(colnames(newZ), expected_Z)

    if (length(missing_Z) > 0L || length(extra_Z) > 0L) {
      stop(
        "The index-component design matrix generated from `newdata` does not match the fitted model.",
        call. = FALSE
      )
    }

    newZ <- newZ[, expected_Z, drop = FALSE]
  }

  ############################################################
  # 12. Subset offset to complete observations
  ############################################################

  offset_complete <- if (is.null(offset)) {
    NULL
  } else {
    as.numeric(offset[complete])
  }

  ############################################################
  # 13. Call existing matrix prediction engine
  ############################################################

  pred_complete <- .predict_pgplsim_matrix(
    object = object,
    newX = newX,
    newZ = newZ,
    type = type,
    offset = offset_complete,
    ...
  )

  ############################################################
  # 14. Restore original newdata row structure
  ############################################################

  out <- rep(NA_real_, nrow(newdata))
  out[complete] <- as.numeric(pred_complete)
  out
}

#' Internal matrix-based prediction engine
#'
#' @keywords internal
#' @noRd
.predict_pgplsim_matrix <- function(
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

  if (!inherits(object, "pgplsim")) {
    stop(
      "object must inherit from class \"pgplsim\".",
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
          "`newX` was not supplied and the fitted object ",
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
          "`newZ` was not supplied and the fitted object ",
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
      "`newX` and `newZ` must have the same number of rows.",
      call. = FALSE
    )
  }

  if (ncol(newX) != length(object$beta_hat)) {
    stop(
      paste0(
        "The linear prediction matrix has ",
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
        "The index prediction matrix has ",
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
        "Prediction covariates must contain only finite, ",
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
            "model using the current pgplsim version."
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

  ############################################################
  # 7. Link- and response-scale prediction
  ############################################################

  family_name <- tolower(
    as.character(
      object$family$family
    )[1L]
  )

  ############################################################
  # Binomial
  ############################################################

  if (identical(family_name, "binomial")) {

    if (!is.null(offset)) {
      stop(
        "`offset` is currently supported only for Poisson prediction.",
        call. = FALSE
      )
    }

    if (identical(type, "link")) {
      return(eta)
    }

    return(
      stats::plogis(eta)
    )
  }

  ############################################################
  # Poisson
  ############################################################

  if (identical(family_name, "poisson")) {

    n_new <- nrow(newX)

    if (is.null(offset)) {

      if (in_sample) {

        if (
          "exposure" %in% names(object) &&
          !is.null(object$exposure)
        ) {

          ######################################################
          # Current public formula fit
          #
          # object$offset:
          #   additive link-scale offset
          #
          # object$exposure:
          #   corresponding positive exposure
          ######################################################

          if (
            "offset" %in% names(object) &&
            !is.null(object$offset)
          ) {

            prediction_offset <- as.numeric(
              object$offset
            )

          } else {

            prediction_offset <- rep(
              0,
              n_new
            )
          }

        } else if (
          "offset" %in% names(object) &&
          !is.null(object$offset)
        ) {

          ######################################################
          # Backward-compatible low-level matrix fit
          #
          # pgplsim_fit() stores positive exposure in
          # object$offset.
          ######################################################

          legacy_exposure <- as.numeric(
            object$offset
          )

          if (
            length(legacy_exposure) != n_new ||
            anyNA(legacy_exposure) ||
            any(!is.finite(legacy_exposure)) ||
            any(legacy_exposure <= 0)
          ) {
            stop(
              "The stored Poisson exposure is invalid.",
              call. = FALSE
            )
          }

          prediction_offset <- log(
            legacy_exposure
          )

        } else {

          ######################################################
          # No fitted offset:
          # exposure = 1, so log-exposure = 0
          ######################################################

          prediction_offset <- rep(
            0,
            n_new
          )
        }

      } else {

        ########################################################
        # New observations without supplied offset:
        # exposure = 1, so log-exposure = 0
        ########################################################

        prediction_offset <- rep(
          0,
          n_new
        )
      }

    } else {

      prediction_offset <- as.numeric(offset)

      if (length(prediction_offset) == 1L) {
        prediction_offset <- rep(
          prediction_offset,
          n_new
        )
      }

      if (length(prediction_offset) != n_new) {
        stop(
          paste0(
            "For Poisson prediction, `offset` must have ",
            "length 1 or the same number of observations ",
            "as the prediction data."
          ),
          call. = FALSE
        )
      }
    }

    if (
      anyNA(prediction_offset) ||
      any(!is.finite(prediction_offset))
    ) {
      stop(
        paste0(
          "For Poisson prediction, `offset` must contain ",
          "only finite, non-missing values."
        ),
        call. = FALSE
      )
    }

    eta_full <- eta + prediction_offset

    if (identical(type, "link")) {
      return(eta_full)
    }

    eta_bounded <- pmin(
      pmax(eta_full, -30),
      30
    )

    return(
      exp(eta_bounded)
    )
  }

  stop(
    paste0(
      "Prediction is not implemented for family: ",
      family_name,
      "."
    ),
    call. = FALSE
  )
}
