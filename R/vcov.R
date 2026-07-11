############################################################
# Variance-covariance methods for LuGPLSIM
############################################################

#' Variance-Covariance Matrix for a LuGPLSIM Model
#'
#' Extracts model-based, sandwich, leverage-adjusted,
#' finite-sample-corrected, or bootstrap variance estimates
#' for a fitted LuGPLSIM model.
#'
#' The internal unconstrained parameter ordering is
#'
#' \deqn{
#'   (\theta^\top, \beta^\top, \gamma^\top)^\top,
#' }
#'
#' where
#'
#' \deqn{
#'   \alpha(\theta)
#'   =
#'   \frac{(1,\theta^\top)^\top}
#'        {\sqrt{1+\|\theta\|^2}}.
#' }
#'
#' The default returned matrix corresponds to the scientifically
#' interpretable parameter vector
#'
#' \deqn{
#'   (\alpha^\top,\beta^\top)^\top.
#' }
#'
#' @param object A fitted object of class `"LuGPLSIM"`.
#' @param type Variance estimator. One of `"model"`, `"sandwich"`,
#'   `"leverage"`, `"finite"`, or `"bootstrap"`.
#' @param component Parameter block to return. One of `"parameters"`,
#'   `"alpha"`, `"beta"`, `"theta"`, `"gamma"`, or `"full"`.
#' @param bootstrap Optional object returned by
#'   [bootstrap_LuGPLSIM()]. Required when `type = "bootstrap"` unless
#'   a bootstrap result is stored in `object$bootstrap`.
#' @param finite_df Optional effective parameter count used for the
#'   finite-sample correction. The default is the dimension of the
#'   unconstrained parameter vector.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A variance-covariance matrix.
#'
#' @export
vcov.LuGPLSIM <- function(
    object,
    type = c(
      "model",
      "sandwich",
      "leverage",
      "finite",
      "bootstrap"
    ),
    component = c(
      "parameters",
      "alpha",
      "beta",
      "theta",
      "gamma",
      "full"
    ),
    bootstrap = NULL,
    finite_df = NULL,
    ...
) {

  if (!inherits(object, "LuGPLSIM")) {
    stop("object must inherit from class 'LuGPLSIM'.", call. = FALSE)
  }

  type <- match.arg(type)
  component <- match.arg(component)

  dims <- lugplsim_parameter_dimensions(object)

  if (type == "bootstrap") {

    if (is.null(bootstrap)) {
      bootstrap <- object$bootstrap
    }

    if (is.null(bootstrap)) {
      stop(
        paste0(
          "A bootstrap result is required. Run\n",
          "  boot <- bootstrap_LuGPLSIM(object, B = 500)\n",
          "and then use\n",
          "  vcov(object, type = 'bootstrap', bootstrap = boot)."
        ),
        call. = FALSE
      )
    }

    return(
      extract_bootstrap_vcov(
        bootstrap = bootstrap,
        component = component,
        object = object
      )
    )
  }

  V_full <- switch(
    type,
    model = lugplsim_model_vcov(object),
    sandwich = lugplsim_sandwich_vcov(
      object,
      adjustment = "none"
    ),
    leverage = lugplsim_sandwich_vcov(
      object,
      adjustment = "leverage"
    ),
    finite = lugplsim_sandwich_vcov(
      object,
      adjustment = "finite",
      finite_df = finite_df
    )
  )

  extract_lugplsim_vcov_component(
    V_full = V_full,
    object = object,
    component = component,
    dims = dims
  )
}


############################################################
# Model-based covariance
############################################################

lugplsim_model_vcov <- function(object) {

  if (!is.null(object$vcov_full)) {
    V <- as.matrix(object$vcov_full)
  } else if (!is.null(object$H_inv)) {
    V <- -as.matrix(object$H_inv)
  } else {
    stop(
      "The fitted object does not contain vcov_full or H_inv.",
      call. = FALSE
    )
  }

  V <- 0.5 * (V + t(V))

  if (any(!is.finite(V))) {
    stop(
      "The model-based covariance matrix contains non-finite values.",
      call. = FALSE
    )
  }

  V
}


############################################################
# Sandwich covariance
############################################################

lugplsim_sandwich_vcov <- function(
    object,
    adjustment = c("none", "leverage", "finite"),
    finite_df = NULL
) {

  adjustment <- match.arg(adjustment)

  bread <- lugplsim_model_vcov(object)
  d <- nrow(bread)

  score <- object$score_contributions

  if (is.null(score)) {
    stop(
      paste0(
        "The fitted object does not contain observation-level ",
        "score contributions.\n",
        "Store an n by d matrix as object$score_contributions ",
        "when fitting the model."
      ),
      call. = FALSE
    )
  }

  score <- as.matrix(score)

  # Permit either n x d or d x n storage.
  if (ncol(score) == d) {
    score_n_by_d <- score
  } else if (nrow(score) == d) {
    score_n_by_d <- t(score)
  } else {
    stop(
      paste0(
        "score_contributions has incompatible dimensions. ",
        "It must be n by ", d, " or ", d, " by n."
      ),
      call. = FALSE
    )
  }

  n <- nrow(score_n_by_d)

  if (any(!is.finite(score_n_by_d))) {
    stop(
      "score_contributions contains non-finite values.",
      call. = FALSE
    )
  }

  if (adjustment == "leverage") {

    h <- lugplsim_hatvalues(object, n = n, d = d)

    # HC2-style adjustment.
    denominator <- sqrt(pmax(1 - h, 1e-8))
    score_n_by_d <- score_n_by_d / denominator
  }

  meat <- crossprod(score_n_by_d)

  V <- bread %*% meat %*% bread
  V <- 0.5 * (V + t(V))

  if (adjustment == "finite") {

    if (is.null(finite_df)) {
      finite_df <- d
    }

    if (!is.numeric(finite_df) ||
        length(finite_df) != 1L ||
        !is.finite(finite_df) ||
        finite_df < 0) {
      stop(
        "finite_df must be one nonnegative finite number.",
        call. = FALSE
      )
    }

    if (n <= finite_df) {
      stop(
        "Finite-sample correction is undefined because n <= finite_df.",
        call. = FALSE
      )
    }

    V <- (n / (n - finite_df)) * V
  }

  V
}


############################################################
# Hat values for leverage-adjusted covariance
############################################################

lugplsim_hatvalues <- function(object, n, d) {

  if (!is.null(object$hatvalues)) {

    h <- as.numeric(object$hatvalues)

    if (length(h) != n) {
      stop(
        "The stored hatvalues vector has the wrong length.",
        call. = FALSE
      )
    }

    return(pmin(pmax(h, 0), 1 - 1e-8))
  }

  E <- object$estimating_design
  w <- object$working_weights

  if (is.null(E) || is.null(w)) {
    stop(
      paste0(
        "Leverage-adjusted covariance requires either object$hatvalues ",
        "or both object$estimating_design and object$working_weights."
      ),
      call. = FALSE
    )
  }

  E <- as.matrix(E)
  w <- as.numeric(w)

  # Expected storage is d x n, matching the likelihood functions.
  if (nrow(E) == d && ncol(E) == n) {
    E_d_by_n <- E
  } else if (ncol(E) == d && nrow(E) == n) {
    E_d_by_n <- t(E)
  } else {
    stop(
      "estimating_design has incompatible dimensions.",
      call. = FALSE
    )
  }

  if (length(w) != n) {
    stop(
      "working_weights has the wrong length.",
      call. = FALSE
    )
  }

  bread <- lugplsim_model_vcov(object)

  EW <- E_d_by_n * rep(sqrt(pmax(w, 0)), each = d)

  # diag(t(EW) %*% bread %*% EW), computed without forming n x n matrix.
  h <- colSums(EW * (bread %*% EW))

  pmin(pmax(as.numeric(h), 0), 1 - 1e-8)
}


############################################################
# Extract covariance blocks
############################################################

extract_lugplsim_vcov_component <- function(
    V_full,
    object,
    component,
    dims = lugplsim_parameter_dimensions(object)
) {

  V_full <- as.matrix(V_full)

  expected_dimension <- dims$p_theta + dims$p_beta + dims$p_gamma

  if (!all(dim(V_full) == expected_dimension)) {
    stop(
      paste0(
        "The full covariance matrix should be ",
        expected_dimension,
        " by ",
        expected_dimension,
        "."
      ),
      call. = FALSE
    )
  }

  theta_index <- seq_len(dims$p_theta)

  beta_index <- dims$p_theta + seq_len(dims$p_beta)

  gamma_index <- dims$p_theta +
    dims$p_beta +
    seq_len(dims$p_gamma)

  if (component == "full") {

    dimnames(V_full) <- list(
      dims$full_names,
      dims$full_names
    )

    return(V_full)
  }

  if (component == "theta") {

    V <- V_full[
      theta_index,
      theta_index,
      drop = FALSE
    ]

    dimnames(V) <- list(
      dims$theta_names,
      dims$theta_names
    )

    return(V)
  }

  if (component == "beta") {

    V <- V_full[
      beta_index,
      beta_index,
      drop = FALSE
    ]

    dimnames(V) <- list(
      dims$beta_names,
      dims$beta_names
    )

    return(V)
  }

  if (component == "gamma") {

    V <- V_full[
      gamma_index,
      gamma_index,
      drop = FALSE
    ]

    dimnames(V) <- list(
      dims$gamma_names,
      dims$gamma_names
    )

    return(V)
  }

  theta_hat <- if (!is.null(object$eta_hat)) {
    as.numeric(object$eta_hat)
  } else {
    theta_from_alpha(object$alpha_hat)
  }

  J_alpha <- jac_alpha_wrt_theta(theta_hat)

  if (component == "alpha") {

    V_theta <- V_full[
      theta_index,
      theta_index,
      drop = FALSE
    ]

    V_alpha <- J_alpha %*% V_theta %*% t(J_alpha)
    V_alpha <- 0.5 * (V_alpha + t(V_alpha))

    dimnames(V_alpha) <- list(
      dims$alpha_names,
      dims$alpha_names
    )

    return(V_alpha)
  }

  # Transform from (theta, beta, gamma) to (alpha, beta).
  T_alpha_beta <- matrix(
    0,
    nrow = dims$p_alpha + dims$p_beta,
    ncol = expected_dimension
  )

  T_alpha_beta[
    seq_len(dims$p_alpha),
    theta_index
  ] <- J_alpha

  T_alpha_beta[
    dims$p_alpha + seq_len(dims$p_beta),
    beta_index
  ] <- diag(dims$p_beta)

  V_parameters <- T_alpha_beta %*%
    V_full %*%
    t(T_alpha_beta)

  V_parameters <- 0.5 * (
    V_parameters + t(V_parameters)
  )

  parameter_names <- c(
    dims$alpha_names,
    dims$beta_names
  )

  dimnames(V_parameters) <- list(
    parameter_names,
    parameter_names
  )

  V_parameters
}


############################################################
# Parameter dimensions and names
############################################################

lugplsim_parameter_dimensions <- function(object) {

  p_alpha <- length(object$alpha_hat)
  p_theta <- p_alpha - 1L
  p_beta <- length(object$beta_hat)
  p_gamma <- length(object$gamma_hat)

  alpha_names <- names(object$alpha_hat)

  if (is.null(alpha_names) ||
      any(!nzchar(alpha_names))) {
    alpha_names <- paste0("alpha", seq_len(p_alpha))
  } else {
    alpha_names <- paste0("alpha_", alpha_names)
  }

  theta_names <- paste0(
    "theta",
    seq_len(p_theta)
  )

  beta_names <- names(object$beta_hat)

  if (is.null(beta_names) ||
      any(!nzchar(beta_names))) {
    beta_names <- paste0("beta", seq_len(p_beta))
  } else {
    beta_names <- paste0("beta_", beta_names)
  }

  gamma_names <- paste0(
    "gamma",
    seq_len(p_gamma)
  )

  list(
    p_alpha = p_alpha,
    p_theta = p_theta,
    p_beta = p_beta,
    p_gamma = p_gamma,
    alpha_names = alpha_names,
    theta_names = theta_names,
    beta_names = beta_names,
    gamma_names = gamma_names,
    full_names = c(
      theta_names,
      beta_names,
      gamma_names
    )
  )
}


############################################################
# Bootstrap covariance extraction
############################################################

extract_bootstrap_vcov <- function(
    bootstrap,
    component,
    object
) {

  if (!inherits(bootstrap, "bootstrap.LuGPLSIM")) {
    stop(
      "bootstrap must inherit from class 'bootstrap.LuGPLSIM'.",
      call. = FALSE
    )
  }

  V <- switch(
    component,
    parameters = bootstrap$vcov_parameters,
    alpha = bootstrap$vcov_alpha,
    beta = bootstrap$vcov_beta,
    theta = bootstrap$vcov_theta,
    gamma = bootstrap$vcov_gamma,
    full = bootstrap$vcov_full
  )

  if (is.null(V)) {
    stop(
      paste0(
        "The requested bootstrap covariance component '",
        component,
        "' is unavailable."
      ),
      call. = FALSE
    )
  }

  as.matrix(V)
}
