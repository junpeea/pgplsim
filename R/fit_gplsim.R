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
#' @param family A GLM family object. Currently supports `binomial()`, `poisson()`,
#'   and `gaussian()`.
#' @param M Number of B-spline basis functions. Default is 80.
#' @param lambda Initial smoothing parameter. Default is 1.
#' @param maxit Maximum number of iterations. Default is 100.
#' @param tol Convergence tolerance. Default is 1e-6.
#'
#' @return An object of class `"LuGPLSIM"`.
#'
#' @importFrom stats binomial gaussian poisson glm.fit
#' @export
LuGPLSIM <- function(
    y,
    X,
    Z,
    family = binomial(),
    M = 80,
    lambda = 1,
    maxit = 100,
    tol = 1e-6
){


  family <- initialize_family(family)


  n <- length(y)
  p <- ncol(Z)


  # initialize alpha

  alpha <- rep(
    1/sqrt(p),
    p
  )


  beta <-
    rep(0,ncol(X))


  for(iter in 1:maxit){


    old <- c(alpha,beta)


    ##################################################
    # single index
    ##################################################

    u <- as.numeric(
      Z %*% alpha
    )


    basis <- compute_basis(
      u,
      M=M
    )


    B <- basis$B
    P <- basis$penalty



    ##################################################
    # penalized GLM update
    ##################################################

    design <-
      cbind(
        X,
        B
      )


    fit <- glm.fit(
      x=design,
      y=y,
      family=family
    )



    coef <- fit$coefficients


    beta <-
      coef[
        seq_len(ncol(X))
      ]


    gamma <-
      coef[
        -(seq_len(ncol(X)))
      ]



    ##################################################
    # update alpha
    # placeholder: your QR update inserted here
    ##################################################

    alpha <- normalize_alpha(alpha)



    ##################################################
    # smoothing
    ##################################################

    lambda <-
      update_lambda_FS(
        gamma,
        P,
        lambda
      )



    if(
      check_convergence(
        old,
        c(alpha,beta),
        tol)
    ) break


  }



  structure(

    list(

      alpha_hat = alpha,

      beta_hat = beta,

      gamma_hat = gamma,

      u_hat = u,

      basis = basis,

      lambda = lambda,

      family = family,

      iterations = iter

    ),

    class="LuGPLSIM"

  )

}
