############################################################
# Nonparametric bootstrap for pgplsim models
############################################################

#' Bootstrap Inference for pgplsim
#'
#' Performs a nonparametric case bootstrap for a fitted `pgplsim`
#' model. For binomial models, a stratified bootstrap can be used
#' to preserve the observed number of observations in each outcome
#' category.
#'
#' Each bootstrap data set is refitted using [pgplsim()].
#'
#' @param object A fitted object of class `"pgplsim"`.
#' @param B Number of successful bootstrap replications. Must be an integer
#'   greater than or equal to 2.
#' @param y Optional original response vector. If omitted,
#'   `object$y` is used.
#' @param X Optional original matrix of linear covariates. If omitted,
#'   `object$X` is used.
#' @param Z Optional original matrix of single-index covariates. If omitted,
#'   `object$Z` is used.
#' @param offset Optional original offset vector. If omitted,
#'   `object$offset` is used.
#' @param method Bootstrap resampling method. `"case"` samples observations
#'   from the full data. `"stratified"` samples separately within outcome
#'   groups and is intended primarily for binomial models.
#' @param seed Optional random-number seed.
#' @param max_attempts Maximum number of attempted fits. Failed bootstrap
#'   fits are replaced until `B` successful fits are obtained or this limit
#'   is reached.
#' @param parallel Logical; use parallel processing when `TRUE`.
#' @param ncores Number of worker processes used when `parallel = TRUE`.
#'   The default is 2. During CRAN checks, the number of workers is capped
#'   at 2.
#' @param trace Logical; print progress information.
#' @param keep_fits Logical; retain successful fitted model objects.
#' @param use_starting_lambda Logical; initialize each bootstrap fit using
#'   the smoothing parameter from the original fitted model.
#' @param ... Additional arguments passed to [pgplsim()]. Arguments supplied
#'   here override corresponding settings recovered from the original fit.
#'
#' @return An object of class `"bootstrap.pgplsim"` containing bootstrap
#'   estimates, covariance matrices, standard errors, and convergence
#'   information.
#'
#' @export
bootstrap_pgplsim <- function(
    object,
    B = 500L,
    y = NULL,
    X = NULL,
    Z = NULL,
    offset = NULL,
    method = c("case", "stratified"),
    seed = NULL,
    max_attempts = max(2L * B, B + 50L),
    parallel = FALSE,
    ncores = 2L,
    trace = interactive(),
    keep_fits = FALSE,
    use_starting_lambda = TRUE,
    ...
) {

  if (!inherits(object, "pgplsim")) {
    stop("object must inherit from class 'pgplsim'.", call. = FALSE)
  }

  method <- match.arg(method)

  B <- validate_bootstrap_integer(
    B,
    name = "B",
    lower = 2L
  )

  max_attempts <- validate_bootstrap_integer(
    max_attempts,
    name = "max_attempts",
    lower = B
  )

  validate_bootstrap_flag(parallel, "parallel")
  validate_bootstrap_flag(trace, "trace")
  validate_bootstrap_flag(keep_fits, "keep_fits")
  validate_bootstrap_flag(use_starting_lambda, "use_starting_lambda")

  if (!is.null(seed)) {
    seed <- validate_bootstrap_integer(
      seed,
      name = "seed",
      lower = 0L
    )
  }

  original_data <- recover_pgplsim_data(
    object = object,
    y = y,
    X = X,
    Z = Z,
    offset = offset
  )

  y <- original_data$y
  X <- original_data$X
  Z <- original_data$Z
  offset <- original_data$offset

  family_name <- object$family$family

  if (method == "stratified" &&
      !identical(family_name, "binomial")) {
    warning(
      paste0(
        "method = 'stratified' is primarily intended for binomial ",
        "models. Proceeding by stratifying on unique response values."
      ),
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Pre-generate resampling indices so serial and parallel runs use
  # the same bootstrap samples.
  index_list <- lapply(
    seq_len(max_attempts),
    function(i) {
      bootstrap_sample_indices(
        y = y,
        method = method
      )
    }
  )

  fit_control <- recover_pgplsim_control(
    object = object,
    use_starting_lambda = use_starting_lambda,
    extra_arguments = list(...)
  )

  dims <- pgplsim_parameter_dimensions(object)

  original_estimates <- extract_pgplsim_estimates(
    object,
    dims = dims
  )

  worker <- function(index) {

    y_b <- y[index]
    X_b <- X[index, , drop = FALSE]
    Z_b <- Z[index, , drop = FALSE]

    offset_b <- if (is.null(offset)) {
      NULL
    } else {
      offset[index]
    }

    args <- c(
      list(
        y = y_b,
        X = X_b,
        Z = Z_b,
        family = object$family,
        offset = offset_b
      ),
      fit_control
    )

    fit_b <- tryCatch(
      do.call(pgplsim_fit, args),
      error = function(e) e
    )

    if (inherits(fit_b, "error")) {

      message(
        "BOOTSTRAP REFIT ERROR: ",
        conditionMessage(fit_b)
      )

      return(
        list(
          success = FALSE,
          converged = FALSE,
          message = conditionMessage(fit_b),
          estimates = NULL,
          fit = NULL
        )
      )
    }

    estimates <- tryCatch(
      extract_pgplsim_estimates(
        fit_b,
        dims = dims
      ),
      error = function(e) NULL
    )

    valid <- !is.null(estimates) &&
      all(
        vapply(
          estimates,
          function(x) all(is.finite(x)),
          logical(1)
        )
      )

    converged <- isTRUE(fit_b$converged)

    list(
      success = valid && converged,
      converged = converged,
      message = if (valid && converged) {
        NA_character_
      } else if (!converged) {
        "Bootstrap fit did not converge."
      } else {
        "Bootstrap fit returned invalid estimates."
      },
      estimates = if (valid) estimates else NULL,
      fit = if (keep_fits && valid) fit_b else NULL
    )
  }

  if (parallel) {

    ncores <- validate_bootstrap_integer(
      ncores,
      name = "ncores",
      lower = 1L
    )

    ncores <- limit_bootstrap_cores(ncores)

    results <- run_pgplsim_bootstrap_parallel(
      index_list = index_list,
      worker = worker,
      ncores = ncores,
      trace = trace
    )

  } else {

    results <- vector("list", max_attempts)
    successful <- 0L

    for (i in seq_len(max_attempts)) {

      results[[i]] <- worker(index_list[[i]])

      if (isTRUE(results[[i]]$success)) {
        successful <- successful + 1L
      }

      if (trace &&
          (i == 1L ||
           i %% 10L == 0L ||
           successful == B)) {
        message(
          "Bootstrap: ",
          successful,
          " successful fits after ",
          i,
          " attempts."
        )
      }

      if (successful >= B) {
        results <- results[seq_len(i)]
        break
      }
    }
  }

  success_indicator <- vapply(
    results,
    function(x) isTRUE(x$success),
    logical(1)
  )

  successful_results <- results[success_indicator]

  if (length(successful_results) < B) {
    stop(
      paste0(
        "Only ",
        length(successful_results),
        " successful bootstrap fits were obtained after ",
        length(results),
        " attempts. Increase max_attempts or inspect convergence."
      ),
      call. = FALSE
    )
  }

  successful_results <- successful_results[seq_len(B)]

  estimate_matrices <- combine_bootstrap_estimates(
    successful_results = successful_results,
    dims = dims
  )

  covariance <- compute_bootstrap_covariances(
    estimate_matrices = estimate_matrices
  )

  standard_errors <- lapply(
    covariance,
    function(V) {
      sqrt(pmax(diag(V), 0))
    }
  )

  failed_messages <- vapply(
    results[!success_indicator],
    function(x) {
      if (is.null(x$message)) {
        NA_character_
      } else {
        x$message
      }
    },
    character(1)
  )

  out <- list(
    call = match.call(),
    original_call = object$call,
    family = object$family,
    method = method,
    B = B,
    attempts = length(results),
    failures = sum(!success_indicator),
    failure_messages = failed_messages,
    original = original_estimates,
    estimates_alpha = estimate_matrices$alpha,
    estimates_beta = estimate_matrices$beta,
    estimates_theta = estimate_matrices$theta,
    estimates_gamma = estimate_matrices$gamma,
    estimates_parameters = estimate_matrices$parameters,
    estimates_full = estimate_matrices$full,
    vcov_alpha = covariance$alpha,
    vcov_beta = covariance$beta,
    vcov_theta = covariance$theta,
    vcov_gamma = covariance$gamma,
    vcov_parameters = covariance$parameters,
    vcov_full = covariance$full,
    se_alpha = standard_errors$alpha,
    se_beta = standard_errors$beta,
    se_theta = standard_errors$theta,
    se_gamma = standard_errors$gamma,
    se_parameters = standard_errors$parameters,
    se_full = standard_errors$full,
    fits = if (keep_fits) {
      lapply(successful_results, `[[`, "fit")
    } else {
      NULL
    },
    seed = seed
  )

  class(out) <- "bootstrap.pgplsim"

  out
}


############################################################
# Internal validation helpers
############################################################

validate_bootstrap_integer <- function(
    x,
    name,
    lower = 0L
) {

  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x != floor(x) ||
      x < lower ||
      x > .Machine$integer.max) {
    stop(
      sprintf(
        "%s must be a single integer greater than or equal to %d.",
        name,
        as.integer(lower)
      ),
      call. = FALSE
    )
  }

  as.integer(x)
}


validate_bootstrap_flag <- function(
    x,
    name
) {

  if (!is.logical(x) ||
      length(x) != 1L ||
      is.na(x)) {
    stop(
      sprintf("%s must be TRUE or FALSE.", name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


limit_bootstrap_cores <- function(ncores) {

  detected <- suppressWarnings(
    parallel::detectCores(logical = FALSE)
  )

  if (length(detected) == 1L &&
      is.finite(detected) &&
      detected >= 1L) {
    ncores <- min(ncores, as.integer(detected))
  }

  # CRAN's incoming checks request that packages limit parallel
  # computation. Respect that request explicitly.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    ncores <- min(ncores, 2L)
  }

  max(1L, ncores)
}


############################################################
# Recover original fitting data
############################################################

recover_pgplsim_data <- function(
    object,
    y = NULL,
    X = NULL,
    Z = NULL,
    offset = NULL
) {

  if (is.null(y)) {
    y <- object$y
  }

  if (is.null(X)) {
    X <- object$X
  }

  if (is.null(Z)) {
    Z <- object$Z
  }

  if (is.null(offset)) {
    offset <- object$offset
  }

  if (is.null(y) || is.null(X) || is.null(Z)) {
    stop(
      paste0(
        "The original y, X, and Z are required for bootstrapping.\n",
        "Either store them in the fitted object or supply them explicitly:\n",
        "  bootstrap_pgplsim(fit, y = y, X = X, Z = Z)"
      ),
      call. = FALSE
    )
  }

  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- length(y)

  if (nrow(X) != n || nrow(Z) != n) {
    stop(
      "The dimensions of y, X, and Z are incompatible.",
      call. = FALSE
    )
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop("y must contain only finite, non-missing values.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("X must contain only finite, non-missing values.", call. = FALSE)
  }

  if (anyNA(Z) || any(!is.finite(Z))) {
    stop("Z must contain only finite, non-missing values.", call. = FALSE)
  }

  if (!is.null(offset)) {

    offset <- as.numeric(offset)

    if (length(offset) != n) {
      stop(
        "offset must have the same length as y.",
        call. = FALSE
      )
    }

    if (anyNA(offset) || any(!is.finite(offset))) {
      stop(
        "offset must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
  }

  list(
    y = y,
    X = X,
    Z = Z,
    offset = offset
  )
}


############################################################
# Recover fitting-control arguments
############################################################

recover_pgplsim_control <- function(
    object,
    use_starting_lambda,
    extra_arguments
) {

  control <- list(
    M = if (!is.null(object$qn)) object$qn else NULL,
    lambda = if (use_starting_lambda &&
                 !is.null(object$lambda)) {
      object$lambda
    } else {
      1
    },
    maxit = if (!is.null(object$control$maxit)) {
      object$control$maxit
    } else {
      100
    },
    tol = if (!is.null(object$control$tol)) {
      object$control$tol
    } else {
      1e-6
    },
    trace = FALSE
  )

  # Remove NULL elements.
  control <- control[
    !vapply(control, is.null, logical(1))
  ]

  # User-supplied arguments override recovered settings.
  if (length(extra_arguments) > 0L) {

    if (is.null(names(extra_arguments)) ||
        any(!nzchar(names(extra_arguments)))) {
      stop(
        "All arguments supplied through ... must be named.",
        call. = FALSE
      )
    }

    control[names(extra_arguments)] <- extra_arguments
  }

  control
}


############################################################
# Generate bootstrap indices
############################################################

bootstrap_sample_indices <- function(
    y,
    method = c("case", "stratified")
) {

  method <- match.arg(method)
  n <- length(y)

  if (method == "case") {
    return(
      sample.int(
        n = n,
        size = n,
        replace = TRUE
      )
    )
  }

  groups <- split(
    seq_len(n),
    factor(y, exclude = NULL)
  )

  sampled <- unlist(
    lapply(
      groups,
      function(index) {
        sample(
          index,
          size = length(index),
          replace = TRUE
        )
      }
    ),
    use.names = FALSE
  )

  # Randomize row ordering after within-group resampling.
  sample(sampled, length(sampled), replace = FALSE)
}


############################################################
# Extract estimates from one fitted model
############################################################

extract_pgplsim_estimates <- function(
    object,
    dims = pgplsim_parameter_dimensions(object)
) {

  alpha <- as.numeric(object$alpha_hat)
  beta <- as.numeric(object$beta_hat)
  gamma <- as.numeric(object$gamma_hat)

  theta <- if (!is.null(object$eta_hat)) {
    as.numeric(object$eta_hat)
  } else {
    theta_from_alpha(alpha)
  }

  if (length(alpha) != dims$p_alpha ||
      length(theta) != dims$p_theta ||
      length(beta) != dims$p_beta ||
      length(gamma) != dims$p_gamma) {
    stop(
      paste0(
        "Bootstrap fit has parameter dimensions different ",
        "from the original fit."
      ),
      call. = FALSE
    )
  }

  names(alpha) <- dims$alpha_names
  names(theta) <- dims$theta_names
  names(beta) <- dims$beta_names
  names(gamma) <- dims$gamma_names

  list(
    alpha = alpha,
    theta = theta,
    beta = beta,
    gamma = gamma,
    parameters = c(alpha, beta),
    full = c(theta, beta, gamma)
  )
}


############################################################
# Combine successful bootstrap estimates
############################################################

combine_bootstrap_estimates <- function(
    successful_results,
    dims
) {

  estimate_names <- c(
    "alpha",
    "theta",
    "beta",
    "gamma",
    "parameters",
    "full"
  )

  out <- lapply(
    estimate_names,
    function(component) {

      matrix_component <- do.call(
        rbind,
        lapply(
          successful_results,
          function(result) {
            result$estimates[[component]]
          }
        )
      )

      matrix_component <- as.matrix(matrix_component)

      if (is.null(colnames(matrix_component))) {
        colnames(matrix_component) <- names(
          successful_results[[1]]$estimates[[component]]
        )
      }

      matrix_component
    }
  )

  names(out) <- estimate_names

  out
}


############################################################
# Compute bootstrap covariance matrices
############################################################

compute_bootstrap_covariances <- function(
    estimate_matrices
) {

  lapply(
    estimate_matrices,
    function(estimates) {

      V <- stats::cov(estimates)
      V <- as.matrix(V)
      V <- 0.5 * (V + t(V))

      V
    }
  )
}


############################################################
# Parallel bootstrap execution
############################################################

run_pgplsim_bootstrap_parallel <- function(
    index_list,
    worker,
    ncores,
    trace
) {

  cluster <- parallel::makeCluster(ncores)

  on.exit(
    parallel::stopCluster(cluster),
    add = TRUE
  )

  # Load the renamed package namespace explicitly on PSOCK workers.
  # loadNamespace() is preferred here to attaching the package with library().
  parallel::clusterCall(
    cluster,
    function(package) {
      loadNamespace(package)
      invisible(NULL)
    },
    "pgplsim"
  )

  results <- parallel::parLapply(
    cluster,
    index_list,
    worker
  )

  if (trace) {

    successful <- sum(
      vapply(
        results,
        function(x) isTRUE(x$success),
        logical(1)
      )
    )

    message(
      "Parallel bootstrap produced ",
      successful,
      " successful fits from ",
      length(results),
      " attempts."
    )
  }

  results
}


############################################################
# Print bootstrap result
############################################################

#' Print a pgplsim Bootstrap Result
#'
#' @param x An object of class `"bootstrap.pgplsim"`.
#' @param digits Number of significant digits.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.bootstrap.pgplsim <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {

  cat("\npgplsim Bootstrap Inference\n")
  cat("===========================\n")
  cat("Method:              ", x$method, "\n", sep = "")
  cat("Successful fits:     ", x$B, "\n", sep = "")
  cat("Attempted fits:      ", x$attempts, "\n", sep = "")
  cat("Failed fits:         ", x$failures, "\n", sep = "")

  cat("\nBootstrap standard errors: alpha\n")
  print(round(x$se_alpha, digits))

  cat("\nBootstrap standard errors: beta\n")
  print(round(x$se_beta, digits))

  invisible(x)
}
