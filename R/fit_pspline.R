############################################################
# Main LuGPLSIM estimator
############################################################

#' Fit a Generalized Partially Linear Single-Index Model
#'
#' Fits a generalized partially linear single-index model of the form
#' g(E[Y | X, Z]) = X beta + phi(Z alpha).
#'
#' @param y Numeric response vector.
#' @param X Numeric matrix of linear covariates.
#' @param Z Numeric matrix of single-index covariates.
#' @param family A GLM family object. Currently supports `binomial()` and
#'   `poisson()`.
#' @param offset Optional offset/exposure vector. For Poisson models, this is
#'   interpreted as an exposure C, so the offset is log(C).
#' @param M Number of inner knots for the spline basis. If `NULL`, uses
#'   `ceiling(n^(1/3))`.
#' @param lambda Initial smoothing parameter. Default is 1.
#' @param maxit Maximum number of Fisher-scoring iterations. Default is 100.
#' @param tol Convergence tolerance. Default is 1e-6.
#' @param trace Logical; if TRUE, prints iteration progress.
#'
#' @return An object of class `"LuGPLSIM"`.
#'
#' @importFrom stats binomial poisson gaussian glm glm.fit plogis quantile
#' @importFrom splines splineDesign
#' @export
LuGPLSIM <- function(
    y,
    X,
    Z,
    family = binomial(),
    offset = NULL,
    M = NULL,
    lambda = 1,
    maxit = 100,
    tol = 1e-6,
    trace = FALSE
) {

  # ------------------------------------------------------------
  # Validate response and design matrices
  # ------------------------------------------------------------

  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- length(y)

  if (nrow(X) != n || nrow(Z) != n) {
    stop(
      "Dimensions of y, X, and Z are incompatible.",
      call. = FALSE
    )
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop(
      "Response y must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  family_name <- family$family

  if (identical(family_name, "binomial")) {
    if (any(!y %in% c(0, 1))) {
      stop(
        "For the binomial family, y must contain only 0 and 1.",
        call. = FALSE
      )
    }
  }

  if (identical(family_name, "poisson")) {
    if (any(y < 0) || any(y != floor(y))) {
      stop(
        paste0(
          "For the Poisson family, y must contain ",
          "non-negative integer counts."
        ),
        call. = FALSE
      )
    }
  }

  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  if (nrow(X) != n || nrow(Z) != n) {
    stop("Dimensions of y, X, and Z are incompatible.")
  }

  family <- initialize_family(family)
  fam <- family$family

  if (!fam %in% c("binomial", "poisson")) {
    stop("The updated Fisher-scoring implementation currently supports binomial() and poisson().")
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

  ############################################################
  # Final likelihood and fitted basis
  ############################################################

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

  final_basis <- build_lugplsim_basis(
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

  if (fam == "binomial") {
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
      family = family,
      offset = if (fam == "poisson") C else NULL,

      # Original data retained for bootstrap refitting
      y = y,
      X = X,
      Z = Z,

      # Fitting controls retained for bootstrap refitting
      control = list(
        M = M,
        qn = qn,
        lambda_start = lambda,
        maxit = maxit,
        tol = tol,
        trace = trace
      ),

      linear_component = linear_component,
      smooth_component = smooth_component,
      fitted_link = linear_predictor,
      linear_predictors = linear_predictor,
      fitted_mean = mu_hat,
      fitted_values = mu_hat,
      residuals = y - mu_hat,

      logLik = last_fit$lik,
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
      working_weights = last_fit$working_weights,
      score_contributions = last_fit$score_contributions,

      iterations = iter,
      converged = converged,
      convergence_difference = theta_diff,
      qn = qn,
      call = match.call(),
      bootstrap = NULL
    ),
    class = "LuGPLSIM"
  )
}


############################################################
# Internal helpers for LuGPLSIM()
############################################################

initialize_beta_gamma <- function(y, X, Z, C, alpha, family, qn) {

  fam <- family$family
  n <- length(y)
  nBeta <- ncol(X)
  nAlpha <- ncol(Z)

  u <- as.numeric(Z %*% alpha)

  if (fam == "binomial") {
    dataset <- cbind(y, X, Z, u)
  } else {
    dataset <- cbind(y, X, Z, C, u)
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

  basis <- build_lugplsim_basis(
    u = u1,
    qn = qn,
    family_name = fam
  )

  BB <- basis$BB

  if (fam == "binomial") {
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


build_lugplsim_basis <- function(u, qn, family_name) {

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
    coeffs <- c(1, rep(-1, m - 1), 1)
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

  basis <- build_lugplsim_basis(U, qn = qn, family_name = "binomial")
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

  basis <- build_lugplsim_basis(U, qn = qn, family_name = "poisson")
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
