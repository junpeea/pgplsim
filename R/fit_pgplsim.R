############################################################
# Internal matrix-based pgplsim fitting engine
############################################################

#' Internal matrix-based fitting engine
#'
#' Fits a generalized partially linear single-index model using
#' preconstructed response and design matrices.
#'
#' This function is called internally by [pgplsim()] and is not part
#' of the public package API.
#'
#' @param y Numeric response vector.
#' @param X Numeric matrix of linear covariates.
#' @param Z Numeric matrix of single-index covariates.
#' @param family A GLM family object.
#' @param offset Optional positive numeric exposure vector for Poisson
#'   models. This is an internal argument: the engine interprets it as
#'   an exposure C and uses `log(C)` as the additive offset. The public
#'   [pgplsim()] interface instead accepts `offset` on the usual
#'   linear-predictor scale.
#' @param M Number of inner knots for the spline basis.
#' @param lambda Initial smoothing parameter.
#' @param maxit Maximum number of Fisher-scoring iterations.
#' @param tol Convergence tolerance.
#' @param trace Logical; if `TRUE`, print iteration progress.
#'
#' @return An object of class `"pgplsim"`.
#'
#' @keywords internal
#' @noRd
pgplsim_fit <- function(
    y,
    X,
    Z,
    family = stats::binomial(),
    offset = NULL,
    M = NULL,
    lambda = 1,
    maxit = 100,
    tol = 1e-6,
    trace = FALSE
){

  # ------------------------------------------------------------
  # Validate response and design matrices
  # ------------------------------------------------------------

  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- length(y)

  if (nrow(X) != n || nrow(Z) != n) {
    stop(
      "Dimensions of `y`, `X`, and `Z` are incompatible.",
      call. = FALSE
    )
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop(
      "`y` must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  family <- initialize_family(family)
  fam <- family$family

  if (identical(fam, "gaussian") && !identical(family$link, "identity")) {
    stop("Gaussian PGPLSIM currently requires the identity link.",
         call. = FALSE)
  }

  if (!fam %in% c("gaussian", "binomial", "poisson")) {
    stop(
      "Only `gaussian()`, `binomial()` and `poisson()` are currently supported.",
      call. = FALSE
    )
  }

  if (fam == "binomial" && any(!y %in% c(0, 1))) {
    stop(
      "For the binomial family, `y` must contain only 0 and 1.",
      call. = FALSE
    )
  }

  if (fam == "poisson" &&
      (any(y < 0) || any(y != floor(y)))) {
    stop(
      paste0(
        "For the Poisson family, `y` must contain ",
        "non-negative integer counts."
      ),
      call. = FALSE
    )
  }

  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  if (fam == "gaussian" && !is.null(offset)) {
    stop("Gaussian offsets are not implemented in Phase 1 Step 1.2.",
         call. = FALSE)
  }

  if (fam == "poisson") {
    if (is.null(offset)) {
      C <- rep(1, n)
    } else {
      C <- as.numeric(offset)
      if (length(C) != n) stop("offset must have the same length as y.")
      if (any(C <= 0)) stop("For Poisson models, offset/exposure values must be positive.")
    }
  } else {
    C <- NULL
  }

  qn <- if (is.null(M)) ceiling(n^(1 / 3)) else M

  ############################################################
  # Step 0: Initial alpha from working GLM
  ############################################################

  if (fam == "gaussian") {
    init_dat <- data.frame(y = y, Z, X)
    init_fit <- stats::lm(y ~ ., data = init_dat)
  }

  if (fam == "binomial") {
    init_dat <- data.frame(y = y, Z, X)
    init_fit <- stats::glm(y ~ ., data = init_dat, family = stats::binomial())
  }

  if (fam == "poisson") {
    init_dat <- data.frame(y = y, Z, X, off = log(C))
    init_fit <- stats::glm(y ~ . - off + offset(off),
                           data = init_dat,
                           family = stats::poisson())
  }

  init_coef <- stats::coef(init_fit)

  alpha <- as.numeric(init_coef[2:(nAlpha + 1)])
  alpha <- alpha / sqrt(sum(alpha^2))

  if (alpha[1] < 0) {
    alpha <- -alpha
  }

  ############################################################
  # Step 1: Initial beta and gamma given initial alpha
  ############################################################

  init_bg <- initialize_beta_gamma(
    y = y,
    X = X,
    Z = Z,
    C = C,
    alpha = alpha,
    family = family,
    qn = qn
  )

  beta <- init_bg$beta
  gamma <- init_bg$gamma

  if (anyNA(gamma)) {
    stop("Initial spline coefficients contain NA values.")
  }

  theta <- c(theta_from_alpha(alpha), beta, gamma)

  lambda_current <- lambda
  converged <- FALSE
  last_fit <- NULL

  ############################################################
  # Fisher-scoring loop
  ############################################################

  # Gaussian: manuscript Section 3.1 uses an INNER parameter fit at fixed
  # smoothing and an OUTER Fellner-Schall update. Keep other families' proven
  # iteration paths unchanged in this phase.
  if (fam == "gaussian") {
    gaussian_stage <- fit_gaussian_two_stage_pkg(
      theta = theta, Y = y, X = X, Z = Z,
      nAlpha = nAlpha, nBeta = nBeta, qn = qn,
      lambda = lambda_current, maxit = maxit, tol = tol, trace = trace
    )
    theta <- gaussian_stage$theta
    alpha <- gaussian_stage$alpha
    beta <- gaussian_stage$beta
    gamma <- gaussian_stage$gamma
    lambda_current <- gaussian_stage$lambda
    converged <- gaussian_stage$converged
    theta_diff <- gaussian_stage$theta_diff
    iter <- gaussian_stage$iterations
    gaussian_smoothing <- gaussian_stage$smoothing
  } else {
    for (iter in seq_len(maxit)) {

    theta_old <- theta
    alpha_old <- alpha
    beta_old <- beta
    gamma_old <- gamma

    if (fam == "binomial") {
      current <- findLogLik_bino_pkg(
        Y = y,
        X = X,
        Z = Z,
        alpha = alpha_old,
        beta = beta_old,
        gamma = gamma_old,
        lambda = lambda_current,
        qn = qn
      )

      theta_new <- FishScore_stephalv_Bino_pkg(
        theta = theta_old,
        lambda = lambda_current,
        mfit = current,
        Y = y,
        X = X,
        Z = Z,
        nAlpha = nAlpha,
        nBeta = nBeta,
        qn = qn
      )
    }

    if (fam == "poisson") {
      current <- findLogLik_pois_pkg(
        Y = y,
        X = X,
        Z = Z,
        C = C,
        alpha = alpha_old,
        beta = beta_old,
        gamma = gamma_old,
        lambda = lambda_current,
        qn = qn
      )

      lambda_current <- update_lambda_pkg(
        mfit = current,
        alpha = alpha_old,
        beta = beta_old,
        gamma = gamma_old,
        lambda = lambda_current
      )

      current <- findLogLik_pois_pkg(
        Y = y,
        X = X,
        Z = Z,
        C = C,
        alpha = alpha_old,
        beta = beta_old,
        gamma = gamma_old,
        lambda = lambda_current,
        qn = qn
      )

      theta_new <- FishScore_stephalv_Pois_pkg(
        theta = theta_old,
        lambda = lambda_current,
        mfit = current,
        Y = y,
        X = X,
        Z = Z,
        C = C,
        nAlpha = nAlpha,
        nBeta = nBeta,
        qn = qn
      )
    }

    if (fam == "binomial") {
      lambda_current <- update_lambda_pkg(
        mfit = current,
        alpha = alpha_old,
        beta = beta_old,
        gamma = gamma_old,
        lambda = lambda_current
      )
    }

    alpha <- alpha_from_theta(theta_new[seq_len(nAlpha - 1)])
    beta <- theta_new[nAlpha:(nAlpha + nBeta - 1)]
    gamma <- theta_new[(nAlpha + nBeta):length(theta_new)]
    theta <- theta_new

    theta_diff <- sqrt(sum((theta - theta_old)^2))

    if (trace) {
      cat("Iteration:", iter,
          " theta_diff:", theta_diff,
          " lambda:", lambda_current, "\n")
    }

    if (theta_diff < tol) {
      converged <- TRUE
      break
    }
  }
  }

  ############################################################
  # Final likelihood and fitted basis
  ############################################################

  if (fam == "gaussian") {
    last_fit <- findLogLik_gaussian_pkg(
      Y = y, X = X, Z = Z, alpha = alpha,
      beta = beta, gamma = gamma,
      lambda = lambda_current, qn = qn
    )
  }

  if (fam == "binomial") {
    last_fit <- findLogLik_bino_pkg(
      Y = y,
      X = X,
      Z = Z,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      lambda = lambda_current,
      qn = qn
    )
  }

  if (fam == "poisson") {
    last_fit <- findLogLik_pois_pkg(
      Y = y,
      X = X,
      Z = Z,
      C = C,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      lambda = lambda_current,
      qn = qn
    )
  }

  final_basis <- build_pgplsim_basis(
    u = as.numeric(Z %*% alpha),
    qn = qn,
    family_name = fam
  )

  u_hat <- final_basis$u
  B_hat <- final_basis$B
  BD_hat <- final_basis$BD
  BB_hat <- final_basis$BB
  BBD_hat <- final_basis$BBD

  knots <- final_basis$knots
  Q2 <- final_basis$Q2
  ord <- final_basis$ord
  outer_ok <- final_basis$outer_ok

  # Save everything required to reproduce the exact transformed
  # spline basis for future observations.
  basis_spec <- list(
    basis_function = "splineDesign",
    knots = knots,
    ord = ord,
    outer_ok = outer_ok,
    Q2 = Q2,
    lower_boundary = min(knots),
    upper_boundary = max(knots),
    n_raw_basis = ncol(B_hat),
    n_transformed_basis = ncol(BB_hat),
    family_name = fam
  )

  if (ncol(BB_hat) != length(gamma)) {
    stop(
      paste0(
        "The final transformed spline basis has ",
        ncol(BB_hat),
        " columns, but gamma has length ",
        length(gamma),
        "."
      ),
      call. = FALSE
    )
  }

  smooth_component <- as.numeric(
    BB_hat %*% gamma
  )

  linear_component <- as.numeric(
    X %*% beta
  )

  linear_predictor <- linear_component +
    smooth_component

  if (fam == "gaussian") {
    mu_hat <- linear_predictor
  } else if (fam == "binomial") {
    mu_hat <- stats::plogis(
      linear_predictor
    )
  } else {
    mu_hat <- C * exp(
      pmin(
        pmax(linear_predictor, -30),
        30
      )
    )
  }

  ############################################################
  # Model-based variance
  ############################################################

  H_inv <- last_fit$H.inv
  V_full <- -H_inv
  if (fam == "gaussian") {
    # H_inv uses the unit-dispersion Fisher curvature; its
    # covariance must be scaled by the estimated residual variance.
    V_full <- last_fit$dispersion * V_full
  }

  p_theta <- nAlpha - 1
  V_theta <- V_full[seq_len(p_theta), seq_len(p_theta), drop = FALSE]

  V_alpha_model <- var_alpha_from_var_theta(
    theta = theta[seq_len(p_theta)],
    V_theta = V_theta
  )

  se_alpha_model <- sqrt(pmax(diag(V_alpha_model), 0))
  se_beta_model <- sqrt(pmax(diag(V_full)[nAlpha:(nAlpha + nBeta - 1)], 0))
  se_eta_model <- sqrt(pmax(diag(V_theta), 0))

  names(alpha) <- colnames(Z)
  names(beta) <- colnames(X)

  structure(
    list(
      alpha_hat = alpha,
      beta_hat = beta,
      gamma_hat = gamma,
      eta_hat = theta[seq_len(p_theta)],
      u_hat = u_hat,

      # Final spline objects corresponding exactly to gamma_hat
      B_hat = B_hat,
      BD_hat = BD_hat,
      BB_hat = BB_hat,
      BBD_hat = BBD_hat,

      # Preferred standardized names for prediction methods
      basis_matrix = BB_hat,
      basis_derivative_matrix = BBD_hat,
      basis_spec = basis_spec,

      # Retained for backward compatibility
      knots = knots,
      Q2 = Q2,
      ord = ord,
      outer_ok = outer_ok,
      lambda = lambda_current,
      smoothing_diagnostics = if (fam == "gaussian") gaussian_smoothing else NULL,
      outer_iterations = if (fam == "gaussian") gaussian_stage$outer_iterations else NULL,
      family = family,
      offset = if (fam == "poisson") C else NULL,

      # Original data retained for bootstrap refitting
      y = y,
      X = X,
      Z = Z,

      linear_component = linear_component,
      smooth_component = smooth_component,

      # Fitting controls retained for bootstrap refitting
      control = list(
        M = M,
        qn = qn,
        lambda_start = lambda,
        maxit = maxit,
        tol = tol,
        trace = trace
      ),

      fitted_link = linear_predictor,
      linear_predictors = linear_predictor,
      fitted_mean = mu_hat,
      fitted_values = mu_hat,
      residuals = y - mu_hat,

      logLik = last_fit$lik,
      # Gaussian diagnostics: `logLik` is evaluated at the adjusted
      # Pearson dispersion, while `logLik_mle` uses RSS / n.  Only the
      # latter maximizes the ordinary Gaussian likelihood conditional
      # on the fitted mean. Neither quantity includes the spline penalty.
      logLik_mle = if (fam == "gaussian") last_fit$lik_mle else NULL,
      dispersion = if (fam == "gaussian") last_fit$dispersion else NULL,
      dispersion_mle = if (fam == "gaussian") last_fit$dispersion_mle else NULL,
      rss = if (fam == "gaussian") last_fit$rss else NULL,
      residual_df = if (fam == "gaussian") last_fit$residual_df else NULL,
      effective_df = if (fam == "gaussian") last_fit$effective_df else NULL,
      penalty_quadratic = if (fam == "gaussian") last_fit$roughness else NULL,
      penalty_value = if (fam == "gaussian") last_fit$penalty else NULL,
      penalized_logLik = last_fit$penalized_lik,
      grad = last_fit$grad,
      H_inv = H_inv,
      vcov_full = V_full,
      vcov_alpha = V_alpha_model,
      vcov_theta = V_theta,
      se_alpha = se_alpha_model,
      se_beta = se_beta_model,
      se_eta = se_eta_model,

      # Observation-level ingredients used by vcov.R
      estimating_design = last_fit$estimating_design,
      working_weights = if (fam == "gaussian") {
        last_fit$working_weights / last_fit$dispersion
      } else {
        last_fit$working_weights
      },
      score_contributions = if (fam == "gaussian") {
        last_fit$score_contributions / last_fit$dispersion
      } else {
        last_fit$score_contributions
      },

      iterations = iter,
      converged = converged,
      convergence_difference = theta_diff,
      qn = qn,
      bootstrap = NULL
    ),
    class = "pgplsim"
  )
}


############################################################
# Internal helpers for pgplsim()
############################################################

initialize_beta_gamma <- function(y, X, Z, C, alpha, family, qn) {

  fam <- family$family
  n <- length(y)
  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  u <- as.numeric(Z %*% alpha)

  if (fam == "poisson") {
    dataset <- cbind(y, X, Z, C, u)
  } else {
    dataset <- cbind(y, X, Z, u)
  }

  dataset <- dataset[order(dataset[, ncol(dataset)]), , drop = FALSE]

  y1 <- dataset[, 1]
  X1 <- dataset[, 2:(nBeta + 1), drop = FALSE]
  Z1 <- dataset[, (nBeta + 2):(nBeta + nAlpha + 1), drop = FALSE]

  if (fam == "poisson") {
    C1 <- dataset[, ncol(dataset) - 1]
  } else {
    C1 <- NULL
  }

  u1 <- dataset[, ncol(dataset)]

  basis <- build_pgplsim_basis(
    u = u1,
    qn = qn,
    family_name = fam
  )

  BB <- basis$BB

  if (fam == "gaussian") {
    fit0 <- stats::lm.fit(x = cbind(X1, BB), y = y1)
  } else if (fam == "binomial") {
    fit0 <- stats::glm.fit(
      x = cbind(X1, BB),
      y = y1,
      family = stats::binomial()
    )
  } else {
    fit0 <- stats::glm.fit(
      x = cbind(X1, BB),
      y = y1,
      family = stats::poisson(),
      offset = log(C1)
    )
  }

  coef0 <- fit0$coefficients

  list(
    beta = coef0[seq_len(ncol(X1))],
    gamma = coef0[-seq_len(ncol(X1))],
    basis = basis
  )
}

#' Control Parameters for pgplsim Fitting
#'
#' Specifies computational control parameters used by [pgplsim()].
#'
#' @param maxit Maximum number of Fisher-scoring iterations.
#'   Default is 100.
#' @param tol Convergence tolerance. Default is `1e-6`.
#' @param trace Logical; if `TRUE`, iteration progress is printed.
#'
#' @return A list containing fitting control parameters.
#'
#' @examples
#' ctrl <- pgplsim_control(
#'   maxit = 200,
#'   tol = 1e-7,
#'   trace = FALSE
#' )
#'
#' ctrl
#'
#' @export
pgplsim_control <- function(
    maxit = 100,
    tol = 1e-6,
    trace = FALSE
) {

  if (length(maxit) != 1L ||
      !is.numeric(maxit) ||
      is.na(maxit) ||
      !is.finite(maxit) ||
      maxit < 1 ||
      maxit != floor(maxit)) {
    stop(
      "`maxit` must be a positive integer.",
      call. = FALSE
    )
  }

  if (length(tol) != 1L ||
      !is.numeric(tol) ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0) {
    stop(
      "`tol` must be a positive finite number.",
      call. = FALSE
    )
  }

  if (length(trace) != 1L ||
      !is.logical(trace) ||
      is.na(trace)) {
    stop(
      "`trace` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  structure(
    list(
      maxit = as.integer(maxit),
      tol = tol,
      trace = trace
    ),
    class = "pgplsim_control"
  )
}

build_pgplsim_basis <- function(u, qn, family_name) {

  u <- as.numeric(u)

  if (family_name == "binomial") {
    trim <- 0.01
    prob <- seq(trim, 1 - trim, length.out = qn)
    psi <- as.vector(unique(stats::quantile(u, prob)))
    boundary_margin <- 0.05 * (max(u) - min(u))
    knots <- c(
      rep(min(u) - boundary_margin, 3),
      psi,
      rep(max(u) + boundary_margin, 3)
    )
    ord <- qn
    outer_ok <- TRUE
  } else {
    prob <- seq(0, 1, length.out = qn)
    psi <- as.vector(unique(stats::quantile(u, prob)))
    knots <- c(rep(min(u), 3), psi, rep(max(u), 3))
    ord <- 4
    outer_ok <- FALSE
  }

  B <- splines::splineDesign(
    knots = knots,
    x = u,
    ord = ord,
    outer.ok = outer_ok
  )

  BD <- splines::splineDesign(
    knots = knots,
    x = u,
    ord = ord,
    derivs = 1,
    outer.ok = outer_ok
  )

  bs <- matrix(colSums(B), nrow = 1)
  res <- qr(t(bs))
  Q <- qr.Q(res, complete = TRUE)
  R <- qr.R(res, complete = TRUE)
  Q2 <- Q[, seq(ncol(R) + 1, ncol(Q)), drop = FALSE]

  BB <- B %*% Q2
  BBD <- BD %*% Q2

  list(
    u = u,
    B = B,
    BD = BD,
    BB = BB,
    BBD = BBD,
    knots = knots,
    Q2 = Q2,
    ord = ord,
    outer_ok = outer_ok
  )
}


generate_P_matrix_pkg <- function(nE, nB, m = 2) {

  num_rows <- nB - m

  D <- matrix(0, nrow = num_rows, ncol = nB)

  for (i in seq_len(num_rows)) {
    indices <- i:(i + m)
    # Order-m finite differences use the signed binomial coefficients.
    # For m = 2 these are (1, -2, 1), annihilating both constant and
    # linear coefficient sequences, as required by the manuscript.
    coeffs <- (-1)^(m - 0:m) * choose(m, 0:m)
    D[i, indices] <- coeffs
  }

  P <- matrix(0, nE, nE)
  P[((nE - nB + 1):nE), ((nE - nB + 1):nE)] <- t(D) %*% D

  P
}


update_lambda_pkg <- function(mfit, alpha, beta, gamma, lambda) {

  theta_n <- c(theta_from_alpha(alpha), beta, gamma)
  nE <- length(theta_n)
  nB <- length(gamma)

  P <- generate_P_matrix_pkg(nE, nB, 2)

  out <- tryCatch(
    {
      Plinv <- MASS::ginv(lambda * P)

      comp1 <- sum(diag(Plinv %*% P))
      comp2 <- sum(diag(-mfit$H.inv %*% P))
      comp3 <- as.numeric(t(theta_n) %*% P %*% theta_n)

      as.numeric((comp1 - comp2) / (comp3 + 1e-8) * lambda)
    },
    error = function(e) {
      warning("lambda was not updated because the penalty matrix was singular.")
      lambda
    }
  )

  if (!is.finite(out) || out <= 0) {
    out <- lambda
  }

  out
}


# Equation (9), expressed in the package's stored penalty multiplier:
#   s = lambda_manuscript^2, P_s = s P
#   s_new = s * [ tr((sP)^- P) - tr((E'E + sP)^-1 P) ] / (theta'Ptheta)
#         = [ rank(P) - s tr((E'E + sP)^-1 P) ] / (theta'Ptheta).
# The rank expression avoids constructing a pseudoinverse of singular sP.
# Both expressions use the unit-dispersion Fisher information; the common
# Gaussian dispersion factor cancels when the penalty is scaled as in Eq (6).
update_lambda_gaussian_FS_pkg <- function(mfit, gamma, lambda,
                                          lower = 1e-8, upper = 1e8) {
  if (!is.finite(lambda) || lambda <= 0) {
    stop("Gaussian smoothing multiplier must be positive and finite.")
  }
  P <- generate_P_matrix_pkg(nrow(mfit$estimating_design), length(gamma), 2)
  Kinv <- -mfit$H.inv
  # A second-difference penalty on gamma has rank length(gamma) - 2.
  rankP <- length(gamma) - 2L
  rough <- as.numeric(crossprod(gamma, P[
    (nrow(P) - length(gamma) + 1L):nrow(P),
    (ncol(P) - length(gamma) + 1L):ncol(P), drop = FALSE
  ] %*% gamma))
  trace_term <- sum(diag(Kinv %*% P))
  numerator <- rankP - lambda * trace_term
  # A zero-roughness component has no informative multiplicative update.
  # Retaining the prior lambda avoids division by zero / arbitrary explosion.
  if (!is.finite(rough) || rough <= .Machine$double.eps ||
      !is.finite(numerator) || numerator <= 0) {
    return(list(lambda = lambda, updated = FALSE,
                numerator = numerator, roughness = rough,
                penalty_rank = rankP, trace_term = trace_term))
  }
  candidate <- numerator / rough
  if (!is.finite(candidate) || candidate <= 0) candidate <- lambda
  candidate <- max(lower, min(upper, candidate))
  list(lambda = candidate, updated = TRUE,
       numerator = numerator, roughness = rough,
       penalty_rank = rankP, trace_term = trace_term)
}

# Two-stage Gaussian algorithm: for each outer iteration, converge the inner
# Fisher-scoring/step-halving updates with smoothing held fixed, then update
# smoothing by the generalized Fellner-Schall rule. `maxit` bounds both loops.
fit_gaussian_two_stage_pkg <- function(theta, Y, X, Z, nAlpha, nBeta,
                                       qn, lambda, maxit, tol, trace) {
  extract <- function(par) {
    list(alpha = alpha_from_theta(par[seq_len(nAlpha - 1L)]),
         beta = par[nAlpha:(nAlpha + nBeta - 1L)],
         gamma = par[(nAlpha + nBeta):length(par)])
  }
  evaluate <- function(par, smoothing) {
    parts <- extract(par)
    findLogLik_gaussian_pkg(
      Y = Y, X = X, Z = Z, alpha = parts$alpha,
      beta = parts$beta, gamma = parts$gamma,
      lambda = smoothing, qn = qn
    )
  }
  theta_diff <- Inf
  relative_lambda_diff <- Inf
  nsteps <- 0L
  converged <- FALSE
  lambda_history <- numeric(0)
  lambda_candidate_history <- numeric(0)
  inner_history <- integer(0)
  fs_updated <- logical(0)
  outer_iterations <- 0L

  for (outer in seq_len(maxit)) {
    outer_iterations <- outer
    lambda_history <- c(lambda_history, lambda)
    theta_before_outer <- theta
    inner_converged <- FALSE
    for (inner in seq_len(maxit)) {
      current <- evaluate(theta, lambda)
      theta_new <- FishScore_stephalv_Gaussian_pkg(
        theta = theta, lambda = lambda, mfit = current,
        Y = Y, X = X, Z = Z,
        nAlpha = nAlpha, nBeta = nBeta, qn = qn
      )
      theta_diff <- sqrt(sum((theta_new - theta)^2))
      theta <- theta_new
      nsteps <- nsteps + 1L
      if (trace) cat("Gaussian outer:", outer, "inner:", inner,
                     "theta_diff:", theta_diff, "lambda:", lambda, "\n")
      if (theta_diff < tol) {
        inner_converged <- TRUE
        break
      }
    }
    inner_history <- c(inner_history, inner)
    fitted <- evaluate(theta, lambda)
    parts <- extract(theta)
    fs <- update_lambda_gaussian_FS_pkg(fitted, parts$gamma, lambda)
    proposed <- fs$lambda
    lambda_candidate_history <- c(lambda_candidate_history, proposed)
    fs_updated <- c(fs_updated, fs$updated)
    relative_lambda_diff <- abs(log(proposed / lambda))
    outer_theta_diff <- sqrt(sum((theta - theta_before_outer)^2))
    if (trace) cat("Gaussian outer:", outer,
                   "lambda proposed:", proposed,
                   "relative log change:", relative_lambda_diff, "\n")
    # Only call the solution converged after the fixed-lambda fit and
    # smoothing update have BOTH stabilized. Keep a lambda at which theta
    # was actually fitted rather than returning an unrefitted proposal.
    if (inner_converged && outer_theta_diff < tol &&
        relative_lambda_diff < tol) {
      converged <- TRUE
      break
    }
    if (outer == maxit) break
    lambda <- proposed
  }
  parts <- extract(theta)
  list(theta = theta, alpha = parts$alpha, beta = parts$beta,
       gamma = parts$gamma, lambda = lambda, iterations = nsteps,
       outer_iterations = outer_iterations, converged = converged,
       theta_diff = theta_diff,
       smoothing = list(lambda_history = lambda_history,
                        proposed_lambda_history = lambda_candidate_history,
                        inner_iterations = inner_history,
                        fs_updated = fs_updated,
                        relative_log_lambda_change = relative_lambda_diff))
}

# Gaussian identity-link specialization of Lu (2026), Section 3.
# Parameter order is (eta, beta, zeta), matching the existing engine.
# The package stores lambda_eff = lambda_manuscript^2.
findLogLik_gaussian_pkg <- function(Y, X, Z, alpha, beta, gamma,
                                    lambda, qn) {
  n <- length(Y)
  nBeta <- ncol(X)
  nAlpha <- ncol(Z)
  U <- as.vector(Z %*% alpha)
  order_u <- order(U)
  Y <- Y[order_u]
  X <- X[order_u, , drop = FALSE]
  Z <- Z[order_u, , drop = FALSE]
  U <- U[order_u]

  basis <- build_pgplsim_basis(U, qn = qn, family_name = "gaussian")
  BB <- basis$BB
  BBD <- basis$BBD
  eta_phi <- as.vector(BBD %*% gamma)
  mu <- as.vector(BB %*% gamma + X %*% beta)
  resid <- Y - mu
  theta_alpha <- theta_from_alpha(alpha)
  k <- length(theta_alpha)
  r2 <- sum(theta_alpha^2)
  norm1 <- 1 / sqrt(1 + r2)
  norm3 <- norm1 / (1 + r2)
  zrest <- Z[, -1, drop = FALSE]
  zlin <- as.vector(Z[, 1] + zrest %*% theta_alpha)
  psi <- norm1 * zrest - norm3 *
    (zlin * matrix(theta_alpha, nrow = n, ncol = k, byrow = TRUE))
  xi <- psi * eta_phi
  E <- rbind(t(xi), t(X), t(BB))
  theta <- c(theta_alpha, beta, gamma)
  P <- generate_P_matrix_pkg(length(theta), ncol(BB), 2)
  Ptheta <- as.vector(P %*% theta)
  roughness <- as.numeric(crossprod(theta, Ptheta))
  rss <- sum(resid^2)
  penalty <- 0.5 * lambda * roughness
  # Score and expected curvature after canceling the common 1/sigma^2.
  grad <- as.vector(E %*% resid - lambda * Ptheta)
  K <- tcrossprod(E) + lambda * P
  K <- 0.5 * (K + t(K))
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) {
    ridge <- 1e-8
    for (j in seq_len(8L)) {
      R <- tryCatch(chol(K + diag(ridge, nrow(K))),
                    error = function(e) NULL)
      if (!is.null(R)) {
        K <- K + diag(ridge, nrow(K))
        break
      }
      ridge <- ridge * 10
    }
  }
  if (is.null(R)) {
    stop("Failed to make the Gaussian Fisher information positive definite.")
  }
  Kinv <- chol2inv(R)
  # Eq. (8) of Lu (2026) calls rho an effective df, but the stated
  # trace is a MODEL effective df.  The adjusted Pearson estimate
  # requires the RESIDUAL df: n - tr[(E'E + lambda P)^-1 E'E].
  # If a numerical ridge is used, K contains that ridge and the trace
  # is computed from the same matrix actually inverted by scoring.
  effective_df <- sum(diag(Kinv %*% tcrossprod(E)))
  residual_df <- n - effective_df
  if (!is.finite(residual_df) || residual_df <= 0) {
    stop("Gaussian residual degrees of freedom must be positive.")
  }
  dispersion <- rss / residual_df
  if (!is.finite(dispersion) || dispersion <= 0) {
    stop("Gaussian residual variance must be positive and finite.")
  }
  lik <- sum(stats::dnorm(Y, mean = mu,
                          sd = sqrt(dispersion), log = TRUE))
  # Distinguish the full density at the adjusted Pearson dispersion
  # from the ordinary Gaussian profile likelihood (conditional on mu),
  # which is maximized at RSS/n.  An ordinary ML profile log likelihood
  # is useful for independent numerical checks, but not an automatic
  # justification for likelihood-ratio tests in this penalized model.
  dispersion_mle <- rss / n
  lik_mle <- sum(stats::dnorm(Y, mean = mu,
                              sd = sqrt(dispersion_mle), log = TRUE))
  # This is the unscaled penalized objective used for scoring/line search.
  objective <- -0.5 * rss - penalty
  list(
    lik = lik,
    lik_mle = lik_mle,
    penalized_lik = lik - penalty / dispersion,
    objective = objective,
    penalty = penalty,
    roughness = roughness,
    rss = rss,
    dispersion_mle = dispersion_mle,
    grad = grad,
    H.inv = -Kinv,
    basis = basis,
    dispersion = dispersion,
    effective_df = effective_df,
    residual_df = residual_df,
    estimating_design = E,
    working_weights = rep(1, n),
    score_contributions = t(E) * resid
  )
}

FishScore_stephalv_Gaussian_pkg <- function(theta, lambda, mfit,
                                            Y, X, Z, nAlpha, nBeta,
                                            qn, maxk = 12L) {
  # Eq. (7): Newton/Fisher step with k = 0, 1, 2, ...;
  # accept only a strict increase of the *penalized* objective.
  direction <- -as.vector(mfit$H.inv %*% mfit$grad)
  for (k in 0:maxk) {
    theta_try <- theta + (0.5)^k * direction
    alpha_try <- alpha_from_theta(theta_try[seq_len(nAlpha - 1L)])
    beta_try <- theta_try[nAlpha:(nAlpha + nBeta - 1L)]
    gamma_try <- theta_try[(nAlpha + nBeta):length(theta_try)]
    candidate <- tryCatch(
      findLogLik_gaussian_pkg(
        Y = Y, X = X, Z = Z, alpha = alpha_try,
        beta = beta_try, gamma = gamma_try,
        lambda = lambda, qn = qn
      ), error = function(e) NULL
    )
    if (!is.null(candidate) && is.finite(candidate$objective) &&
        candidate$objective > mfit$objective) {
      return(theta_try)
    }
  }
  # A stationary point or unsuccessful line search leaves theta unchanged.
  theta
}


FishScore_stephalv_Bino_pkg <- function(
    theta,
    lambda,
    mfit,
    Y,
    X,
    Z,
    nAlpha,
    nBeta,
    qn,
    maxk = 8
) {

  grad <- mfit$grad
  H_inv <- mfit$H.inv
  lik_diff <- -Inf
  theta_new <- theta

  for (k in seq_len(maxk)) {

    theta_try <- theta - (0.5)^k * as.vector(H_inv %*% grad)

    alpha_try <- alpha_from_theta(theta_try[seq_len(nAlpha - 1)])
    beta_try <- theta_try[nAlpha:(nAlpha + nBeta - 1)]
    gamma_try <- theta_try[(nAlpha + nBeta):length(theta_try)]

    new <- findLogLik_bino_pkg(
      Y = Y,
      X = X,
      Z = Z,
      alpha = alpha_try,
      beta = beta_try,
      gamma = gamma_try,
      lambda = lambda,
      qn = qn
    )

    lik_diff <- new$lik - mfit$lik

    if (lik_diff > 0) {
      theta_new <- theta_try
      break
    }
  }

  theta_new
}


FishScore_stephalv_Pois_pkg <- function(
    theta,
    lambda,
    mfit,
    Y,
    X,
    Z,
    C,
    nAlpha,
    nBeta,
    qn,
    maxk = 8
) {

  grad <- mfit$grad
  H_inv <- mfit$H.inv
  lik_diff <- -Inf
  theta_new <- theta

  for (k in seq_len(maxk)) {

    theta_try <- theta - (0.5)^k * as.vector(H_inv %*% grad)

    alpha_try <- alpha_from_theta(theta_try[seq_len(nAlpha - 1)])
    beta_try <- theta_try[nAlpha:(nAlpha + nBeta - 1)]
    gamma_try <- theta_try[(nAlpha + nBeta):length(theta_try)]

    new <- findLogLik_pois_pkg(
      Y = Y,
      X = X,
      Z = Z,
      C = C,
      alpha = alpha_try,
      beta = beta_try,
      gamma = gamma_try,
      lambda = lambda,
      qn = qn
    )

    lik_diff <- new$lik - mfit$lik

    if (lik_diff > 0) {
      theta_new <- theta_try
      break
    }
  }

  theta_new
}


findLogLik_bino_pkg <- function(Y, X, Z, alpha, beta, gamma, lambda, qn) {

  n <- length(Y)
  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  U <- as.vector(Z %*% alpha)
  dataset <- cbind(Y, X, Z, U)
  dataset <- dataset[order(dataset[, ncol(dataset)]), , drop = FALSE]

  Y <- dataset[, 1]
  X <- dataset[, 2:(nBeta + 1), drop = FALSE]
  Z <- dataset[, (nBeta + 2):(nBeta + nAlpha + 1), drop = FALSE]
  U <- dataset[, ncol(dataset)]

  basis <- build_pgplsim_basis(U, qn = qn, family_name = "binomial")
  BB <- basis$BB
  BBD <- basis$BBD

  eta_phi <- as.vector(BBD %*% gamma)

  lin_pred <- as.vector(BB %*% gamma + X %*% beta)
  mu <- stats::plogis(lin_pred)
  mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

  theta_alpha <- theta_from_alpha(alpha)
  k <- length(theta_alpha)
  r2 <- sum(theta_alpha^2)

  norm1 <- 1 / sqrt(1 + r2)
  norm3 <- norm1 / (1 + r2)

  z1 <- Z[, 1]
  zrest <- Z[, -1, drop = FALSE]
  zlin <- as.vector(z1 + zrest %*% theta_alpha)

  psi_mat <- norm1 * zrest -
    norm3 * (zlin * matrix(theta_alpha,
                           nrow = nrow(Z),
                           ncol = k,
                           byrow = TRUE))

  E <- rbind(
    t(psi_mat) %*% diag(eta_phi),
    t(X),
    t(BB)
  )

  P <- generate_P_matrix_pkg(nrow(E), ncol(BB), 2)

  theta <- c(theta_alpha, beta, gamma)

  lik <- sum(Y * log(mu) + (1 - Y) * log(1 - mu))

  penalty <- 0.5 * lambda *
    as.numeric(crossprod(theta, P %*% theta))
  penalized_lik <- lik - penalty

  residual <- Y - mu

  # One row per observation; columns follow (theta, beta, gamma).
  # The likelihood routine sorts observations by the current index, but
  # cross-products used by sandwich estimators are invariant to row order.
  score_contributions <- t(E) * residual

  grad <- E %*% residual - lambda * P %*% theta

  w <- pmax(mu * (1 - mu), 1e-8)
  Es <- E * rep(sqrt(w), each = nrow(E))

  K <- tcrossprod(Es) + lambda * P
  K <- 0.5 * (K + t(K))

  R <- tryCatch(chol(K), error = function(e) NULL)

  if (is.null(R)) {
    ridge <- 1e-8
    for (j in seq_len(8)) {
      R <- tryCatch(chol(K + diag(ridge, nrow(K))),
                    error = function(e) NULL)
      if (!is.null(R)) {
        K <- K + diag(ridge, nrow(K))
        break
      }
      ridge <- ridge * 10
    }
  }

  if (is.null(R)) {
    stop("Failed to make the binomial Hessian positive definite.")
  }

  Kinv <- chol2inv(R)
  H_inv <- -Kinv

  list(
    lik = lik,
    penalized_lik = penalized_lik,
    penalty = penalty,
    grad = grad,
    H.inv = H_inv,
    basis = basis,
    estimating_design = E,
    working_weights = w,
    score_contributions = score_contributions
  )
}


findLogLik_pois_pkg <- function(Y, X, Z, C, alpha, beta, gamma, lambda, qn) {

  n <- length(Y)
  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  U <- as.vector(Z %*% alpha)
  dataset <- cbind(Y, X, Z, C, U)
  dataset <- dataset[order(dataset[, ncol(dataset)]), , drop = FALSE]

  Y <- dataset[, 1]
  X <- dataset[, 2:(nBeta + 1), drop = FALSE]
  Z <- dataset[, (nBeta + 2):(nBeta + nAlpha + 1), drop = FALSE]
  C <- dataset[, ncol(dataset) - 1]
  U <- dataset[, ncol(dataset)]

  basis <- build_pgplsim_basis(U, qn = qn, family_name = "poisson")
  BB <- basis$BB
  BBD <- basis$BBD

  eta_phi <- as.vector(BBD %*% gamma)

  eta <- as.vector(BB %*% gamma + X %*% beta)
  eta <- pmin(pmax(eta, -30), 30)

  mu <- C * exp(eta)
  mu <- pmax(mu, 1e-8)

  theta_alpha <- theta_from_alpha(alpha)
  k <- length(theta_alpha)
  r2 <- sum(theta_alpha^2)

  norm1 <- 1 / sqrt(1 + r2)
  norm3 <- norm1 / (1 + r2)

  z1 <- Z[, 1]
  zrest <- Z[, -1, drop = FALSE]
  zlin <- as.vector(z1 + zrest %*% theta_alpha)

  psi <- norm1 * zrest -
    norm3 * (zlin * matrix(theta_alpha,
                           nrow = nrow(Z),
                           ncol = k,
                           byrow = TRUE))

  delta <- Y - mu

  E <- rbind(
    t(psi) %*% diag(eta_phi),
    t(X),
    t(BB)
  )

  P <- generate_P_matrix_pkg(nrow(E), ncol(BB), 2)

  theta <- c(theta_alpha, beta, gamma)

  lik <- sum(Y * log(mu) - mu)

  penalty <- 0.5 * lambda *
    as.numeric(crossprod(theta, P %*% theta))
  penalized_lik <- lik - penalty

  # One row per observation; columns follow (theta, beta, gamma).
  score_contributions <- t(E) * delta

  grad <- E %*% delta - lambda * P %*% theta

  Es <- E * rep(sqrt(mu), each = nrow(E))
  K <- tcrossprod(Es) + lambda * P
  K <- 0.5 * (K + t(K))

  R <- tryCatch(chol(K), error = function(e) NULL)

  if (is.null(R)) {
    ridge <- 1e-8
    for (j in seq_len(8)) {
      R <- tryCatch(chol(K + diag(ridge, nrow(K))),
                    error = function(e) NULL)
      if (!is.null(R)) {
        K <- K + diag(ridge, nrow(K))
        break
      }
      ridge <- ridge * 10
    }
  }

  if (is.null(R)) {
    stop("Failed to make the Poisson Hessian positive definite.")
  }

  Kinv <- chol2inv(R)
  H_inv <- -Kinv

  list(
    lik = lik,
    penalized_lik = penalized_lik,
    penalty = penalty,
    grad = grad,
    H.inv = H_inv,
    basis = basis,
    estimating_design = E,
    working_weights = mu,
    score_contributions = score_contributions
  )
}

#' Fit a generalized partially linear single-index model
#'
#' Fits a generalized partially linear single-index model of the form
#'
#' \deqn{
#' g\{E(Y_i \mid X_i, Z_i)\}
#' =
#' X_i^\top \beta
#' +
#' \phi(Z_i^\top \alpha),
#' }
#'
#' where \eqn{g(\cdot)} is the link function, \eqn{X_i} contains
#' covariates with linear effects, \eqn{Z_i} contains covariates entering
#' the nonlinear single-index component, \eqn{\beta} is the vector of
#' linear coefficients, \eqn{\alpha} is the single-index coefficient
#' vector, and \eqn{\phi(\cdot)} is an unknown smooth function estimated
#' using spline basis functions.
#'
#' The single-index coefficient vector is normalized to have unit Euclidean
#' norm. Its sign is oriented so that the first component is positive,
#' providing an identifiable representation of the index.
#'
#' @param formula A two-sided model formula specifying the response and the
#'   linear component of the model. A standard intercept is included unless
#'   explicitly removed in the formula. For example, `y ~ x1 + x2`.
#'
#' @param index A one-sided formula specifying the numeric covariates entering
#'   the nonlinear single-index component. For example, `~ z1 + z2 + z3`.
#'   An intercept is not included in the index.
#'
#' @param data A data frame containing the variables referenced in `formula`
#'   and `index`.
#'
#' @param family A family object describing the response distribution and
#'   link function. Currently, `binomial()` and `poisson()` are supported.
#'
#' @param offset Optional numeric offset on the linear-predictor scale for
#'   Poisson models. For an exposure such as population or person-time,
#'   supply its logarithm; for example, `offset = log(population)`.
#'   The offset must have the same length as the number of rows in `data`.
#'   Missing offset values are handled jointly with missing values in the
#'   model variables. Offsets are currently not supported for binomial models.
#'
#' @param M Optional positive integer controlling the spline basis
#'   construction. If `NULL`, the default is
#'   `ceiling(n^(1 / 3))`, where `n` is the analysis-sample size.
#'
#' @param lambda Positive initial smoothing parameter controlling the penalty
#'   applied to the nonlinear component. The smoothing parameter may be
#'   updated during model fitting. Default is `1`.
#'
#' @param control A control object created by [pgplsim_control()] specifying
#'   convergence and tracing options.
#'
#' @details
#' The model combines a parametric linear component with a flexible nonlinear
#' function of a single index. The linear design matrix is constructed from
#' `formula` using standard R model-formula conventions, including factor
#' coding and contrasts. The single-index design matrix is constructed from
#' `index`; index variables must currently be numeric.
#'
#' Rows containing missing values in the response, linear predictors,
#' single-index predictors, or supplied offset are omitted using one common
#' analysis sample.
#'
#' For Poisson models with an offset \eqn{o_i}, the fitted model is
#'
#' \deqn{
#' \log\{E(Y_i \mid X_i, Z_i)\}
#' =
#' o_i
#' +
#' X_i^\top \beta
#' +
#' \phi(Z_i^\top \alpha).
#' }
#'
#' Thus, when \eqn{E_i} is a population or person-time exposure, specify
#' `offset = log(E_i)`.
#'
#' @return An object of class `"pgplsim"`. The returned object contains,
#'   among other components:
#'
#' \describe{
#'   \item{alpha_hat}{Estimated normalized single-index coefficients.}
#'   \item{beta_hat}{Estimated coefficients for the linear component.}
#'   \item{gamma_hat}{Estimated spline coefficients for the nonlinear
#'     single-index component.}
#'   \item{eta_hat}{Estimated unconstrained parameter vector used to
#'     parameterize the normalized single-index coefficient vector
#'     \eqn{\alpha}.}
#'   \item{u_hat}{Estimated single-index values
#'     \eqn{Z_i^\top \hat{\alpha}}.}
#'   \item{linear_component}{Estimated linear contribution
#'     \eqn{X_i^\top \hat{\beta}}.}
#'   \item{smooth_component}{Estimated nonlinear contribution
#'     \eqn{\hat{\phi}(Z_i^\top \hat{\alpha})}.}
#'   \item{fitted_link}{Estimated observation-level linear predictor,
#'     excluding any Poisson offset.}
#'   \item{fitted_mean}{Estimated response means. For Poisson models,
#'     the supplied offset is incorporated in these fitted means.}
#'   \item{logLik}{Model log-likelihood evaluated at the fitted parameters.}
#'   \item{penalized_logLik}{Penalized model log-likelihood.}
#'   \item{converged}{Logical indicator of algorithmic convergence.}
#'   \item{iterations}{Number of fitting iterations used.}
#'   \item{lambda}{Final smoothing parameter used by the fitted model.}
#'   \item{formula}{The supplied linear-component formula.}
#'   \item{index}{The supplied single-index formula.}
#'   \item{X}{Linear-component model matrix used for fitting.}
#'   \item{Z}{Single-index model matrix used for fitting.}
#'   \item{y}{Response vector used for fitting.}
#'   \item{offset}{For Poisson models, the link-scale offset used for
#'     fitting, or `NULL` if none was supplied.}
#'   \item{exposure}{For Poisson models, `exp(offset)` when an offset was
#'     supplied, and a vector of ones otherwise.}
#'   \item{control}{A list containing the fitting settings used by the
#'     numerical algorithm, including spline-basis and convergence controls.}
#' }
#'
#' @seealso
#' [pgplsim_control()], [predict.pgplsim()], [summary.pgplsim()],
#' [vcov.pgplsim()], [bootstrap_pgplsim()]
#'
#' @examples
#' set.seed(2026)
#'
#' n <- 150
#'
#' dat <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   z1 = rnorm(n),
#'   z2 = rnorm(n)
#' )
#'
#' eta <- 0.5 * dat$x1 -
#'   0.3 * dat$x2 +
#'   sin(dat$z1 + dat$z2)
#'
#' dat$y <- rbinom(
#'   n,
#'   size = 1,
#'   prob = plogis(eta)
#' )
#'
#' fit <- pgplsim(
#'   y ~ x1 + x2,
#'   index = ~ z1 + z2,
#'   data = dat,
#'   family = stats::binomial(),
#'   M = 4
#' )
#'
#' summary(fit)
#'
#' predict(
#'   fit,
#'   newdata = dat[1:5, ],
#'   type = "response"
#' )
#'
#' @export
pgplsim <- function(
    formula,
    index,
    data,
    family = stats::binomial(),
    offset = NULL,
    M = NULL,
    lambda = 1,
    control = pgplsim_control()
) {

  cl <- match.call()

  ############################################################
  # Validate public interface
  ############################################################

  if (!inherits(formula, "formula")) {
    stop(
      "`formula` must be a model formula.",
      call. = FALSE
    )
  }

  if (length(formula) != 3L) {
    stop(
      "`formula` must be a two-sided formula such as `y ~ x1 + x2`.",
      call. = FALSE
    )
  }

  if (!inherits(index, "formula")) {
    stop(
      "`index` must be a one-sided formula.",
      call. = FALSE
    )
  }

  if (length(index) != 2L) {
    stop(
      "`index` must be a one-sided formula such as `~ z1 + z2`.",
      call. = FALSE
    )
  }

  if (!is.data.frame(data)) {
    stop(
      "`data` must be a data frame.",
      call. = FALSE
    )
  }

  if (!inherits(control, "pgplsim_control")) {
    stop(
      "`control` must be created by `pgplsim_control()`.",
      call. = FALSE
    )
  }

  ############################################################
  # Initialize and validate family
  ############################################################

  family <- initialize_family(family)
  fam <- family$family

  if (identical(fam, "gaussian") && !identical(family$link, "identity")) {
    stop("Gaussian PGPLSIM currently requires the identity link.",
         call. = FALSE)
  }

  if (!fam %in% c("gaussian", "binomial", "poisson")) {
    stop(
      "Only `gaussian()`, `binomial()` and `poisson()` are currently supported.",
      call. = FALSE
    )
  }

  ############################################################
  # Construct model frames using standard R formula semantics
  ############################################################

  formula_terms <- stats::terms(
    formula,
    data = data
  )

  index_terms <- stats::terms(
    index,
    data = data
  )

  if (attr(formula_terms, "response") != 1L) {
    stop(
      "`formula` must contain exactly one response.",
      call. = FALSE
    )
  }

  if (attr(index_terms, "response") != 0L) {
    stop(
      "`index` must be a one-sided formula.",
      call. = FALSE
    )
  }

  mf_linear_full <- stats::model.frame(
    formula = formula_terms,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )

  mf_index_full <- stats::model.frame(
    formula = index_terms,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )

  if (nrow(mf_linear_full) != nrow(data) ||
      nrow(mf_index_full) != nrow(data)) {
    stop(
      paste0(
        "Model-frame construction produced incompatible ",
        "numbers of observations."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # Validate index variables
  ############################################################

  index_classes <- vapply(
    mf_index_full,
    function(x) {
      is.numeric(x) || is.integer(x)
    },
    logical(1)
  )

  if (!all(index_classes)) {
    bad_index_vars <- names(mf_index_full)[!index_classes]

    stop(
      paste0(
        "Variables in `index` must currently be numeric. ",
        "Non-numeric variable(s): ",
        paste(bad_index_vars, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  ############################################################
  # Determine one common analysis sample
  ############################################################

  complete <- stats::complete.cases(mf_linear_full) &
    stats::complete.cases(mf_index_full)

  ############################################################
  # Validate public offset
  ############################################################

  if (fam != "poisson" && !is.null(offset)) {
    stop(
      "`offset` is currently supported only for Poisson models.",
      call. = FALSE
    )
  }

  if (!is.null(offset)) {

    if (!is.numeric(offset)) {
      stop(
        "`offset` must be a numeric vector.",
        call. = FALSE
      )
    }

    if (length(offset) != nrow(data)) {
      stop(
        paste0(
          "`offset` must have the same length as the number ",
          "of rows in `data`."
        ),
        call. = FALSE
      )
    }

    if (any(is.infinite(offset))) {
      stop(
        "`offset` must not contain infinite values.",
        call. = FALSE
      )
    }

    complete <- complete & !is.na(offset)
  }

  if (!any(complete)) {
    stop(
      "No complete observations remain after removing missing values.",
      call. = FALSE
    )
  }

  ############################################################
  # Restrict model frames to identical observations
  ############################################################

  mf_linear <- mf_linear_full[
    complete,
    ,
    drop = FALSE
  ]

  mf_index <- mf_index_full[
    complete,
    ,
    drop = FALSE
  ]

  ############################################################
  # Response
  ############################################################

  y <- stats::model.response(mf_linear)

  ############################################################
  # Linear design matrix
  ############################################################

  X <- stats::model.matrix(
    object = stats::delete.response(formula_terms),
    data = mf_linear
  )

  ############################################################
  # Single-index design matrix
  #
  # The single-index component never contains an intercept.
  ############################################################

  index_terms_no_intercept <- stats::terms(
    stats::update.formula(
      index,
      ~ . - 1
    ),
    data = data
  )

  Z <- stats::model.matrix(
    object = index_terms_no_intercept,
    data = mf_index
  )

  if (ncol(Z) < 2L) {
    stop(
      "`index` must contain at least two numeric predictors.",
      call. = FALSE
    )
  }

  ############################################################
  # Restrict public offset to the analysis sample
  ############################################################

  if (!is.null(offset)) {
    offset <- as.numeric(offset[complete])
  }

  ############################################################
  # Convert public link-scale offset to internal exposure scale
  ############################################################

  internal_offset <- if (fam == "poisson") {

    if (is.null(offset)) {

      NULL

    } else {

      C <- exp(offset)

      if (any(!is.finite(C))) {
        stop(
          paste0(
            "`offset` contains values that are too large to ",
            "convert safely to the exposure scale."
          ),
          call. = FALSE
        )
      }

      if (any(C <= 0)) {
        stop(
          paste0(
            "`offset` contains values that are too small to ",
            "convert safely to a positive exposure."
          ),
          call. = FALSE
        )
      }

      C
    }

  } else {
    NULL
  }

  ############################################################
  # Fit using internal matrix engine
  ############################################################

  fit <- pgplsim_fit(
    y = y,
    X = X,
    Z = Z,
    family = family,
    offset = internal_offset,
    M = M,
    lambda = lambda,
    maxit = control$maxit,
    tol = control$tol,
    trace = control$trace
  )

  ############################################################
  # Store public offset semantics
  ############################################################

  if (fam == "poisson") {

    fit["offset"] <- list(offset)

    fit["exposure"] <- list(
      if (is.null(internal_offset)) {
        rep(1, length(y))
      } else {
        internal_offset
      }
    )

  } else {

    # Use list assignment so these components are retained
    # explicitly as NULL rather than deleted from the fitted object.
    fit["offset"] <- list(NULL)
    fit["exposure"] <- list(NULL)
  }

  fit$offset_supplied <- !is.null(offset)

  ############################################################
  # Add formula-interface metadata
  ############################################################

  fit$call <- cl
  fit$formula <- formula
  fit$index <- index

  fit$terms <- formula_terms
  fit$index_terms <- index_terms_no_intercept

  fit$contrasts <- attr(X, "contrasts")

  fit$xlevels <- stats::.getXlevels(
    formula_terms,
    mf_linear
  )

  fit$index_contrasts <- attr(Z, "contrasts")

  fit$index_xlevels <- stats::.getXlevels(
    index_terms_no_intercept,
    mf_index
  )

  fit$nobs <- sum(complete)

  fit$na.action <- if (all(complete)) {
    NULL
  } else {
    structure(
      which(!complete),
      class = "omit"
    )
  }

  fit$data_names <- names(data)

  fit
}
