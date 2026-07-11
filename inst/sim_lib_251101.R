library(matrixcalc)
library(splines)
library(MASS)
library(psych)
library(lqmm)

compute_bmatrix = function(u,B,n){
  
  # Create adjusted design matrix (proposed)
  # Normalize rows of the design matrix using QR decomposition to accommodate sum-to-one constraints
  bs      <- matrix(colSums(B),nrow=1)
  res     <- qr(t(bs))
  bs.Q    <- qr.Q(res,complete = TRUE)
  bs.R    <- qr.R(res,complete = TRUE)
  bs.Q2   <- bs.Q[,seq(ncol(bs.R)+1,ncol(bs.Q))]
  matrix  <- B%*%bs.Q2   ## new design matrix compared to original B
  
  return(matrix)
}
# compute_bmatrix(u.1,B.1,n)

generate_P_matrix <- function(nE, nB, m=2) {
  # Number of rows in D
  num_rows <- nB - m
  
  # Initialize an empty matrix
  D <- matrix(0, nrow = num_rows, ncol = nB)
  
  # Fill in the difference operator rows
  for (i in 1:num_rows) {
    # Define the indices for the difference coefficients
    indices <- i:(i + m)
    # Define the coefficients for the difference operator
    coeffs <- c(1, rep(-1, m - 1), 1) # Example: m = 2 -> [1, -2, 1]
    # Assign coefficients to the appropriate row
    D[i, indices] <- coeffs
  }
  
  P = matrix(0,nE,nE)
  # P[(1:nB),(1:nB)] <- t(D) %*% D
  P[((nE-nB+1):nE),((nE-nB+1):nE)] <- t(D) %*% D
  
  return(P)
}

update_lambda = function(mfit,alpha,beta,gamma,lambda){
  
  theta.n<-c(theta_from_alpha(alpha),beta,gamma)
  nE = length(theta.n)
  nB = length(gamma)
  P = generate_P_matrix(nE,nB,2)
  Plinv = ginv(lambda*P)
  
  tryCatch(
    {
      comp1 = tr(Plinv%*%P)
      comp2 = tr(-mfit$H.inv%*%P)
      comp3 = t(theta.n) %*% P %*% theta.n
      return(as.numeric((comp1-comp2)/comp3*lambda))
    },
    error = function(e) {
      message("Matrix is computationally singular. lambda has not been updated at this time")
      return(lambda)
    }
  )
}

FishScore_stephalv_Bino = function(theta,lambda,mfit,maxk=8){
  
  grad<-mfit$grad
  H.inv<-mfit$H.inv
  for(k in 1:maxk){
    theta.n<-theta-(0.5)^k*H.inv%*%grad  
    # alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    # if(sum(is.na(alpha.n))){
    #   alpha.n <- c(0,theta.n[1:(nAlpha-1)])
    #   break;
    # }
    alpha.n<-alpha_from_theta(theta.n[1:(nAlpha-1)]) # version 2025
    beta.n<-theta.n[nAlpha:(nAlpha+nBeta-1)]
    gamma.n<-theta.n[(nBeta+nAlpha):length(theta.n)]
    new<-findLogLik_bino(Y,X,Z,alpha.n,beta.n,gamma.n,lambda)
    lik.diff<-new$lik-mfit$lik
    if(lik.diff>0) break;
  }
  if(lik.diff<=0) theta.n <- theta
  return(theta.n)
  
}

"findLogLik_bino"<-function(Y,X,Z,alpha,beta,gamma,lambda){
  n = length(Y)
  U<-as.vector(Z%*%alpha)
  dataset<-cbind(Y,X,Z,U)
  dataset<-dataset[order(dataset[,ncol(dataset)]),]
  nBeta = ncol(X); nAlpha = ncol(Z)
  Y<-dataset[,1];X<-dataset[,2:(nBeta+1)];Z<-dataset[,(nBeta+2):(nBeta+nAlpha+1)]
  U<-dataset[,ncol(dataset)]
  
  #-----------Find B and derivative of B----------------  
  qn    <- ceiling(n^(1/3))                 ## No. of inner knots
  # qn    <- 4
  trim <- 0.01  # remove extreme values
  prob <- seq(trim, 1 - trim, length.out = qn)
  psi   <- as.vector(unique(quantile(U,prob)))
  boundary_margin <- 0.05 * (max(U) - min(U))
  boundary_knots <- c(rep((min(U) - boundary_margin),3),psi,rep(max(U) + boundary_margin,3))
  B    <- splineDesign(knots = boundary_knots,U,ord = qn,outer.ok = TRUE) 
  BD   <- splineDesign(knots = boundary_knots,U,ord = qn,outer.ok = TRUE,derivs=1)
  
  bs      <- matrix(colSums(B),nrow=1)
  res     <- qr(t(bs))
  bs.Q    <- qr.Q(res,complete = TRUE)
  bs.R    <- qr.R(res,complete = TRUE)
  bs.Q2   <- bs.Q[,seq(ncol(bs.R)+1,ncol(bs.Q))]
  BB   <- B%*%bs.Q2   ## new design matrix compared to original B 
  BBD  <- BD%*%bs.Q2   ## new design matrix compared to original B 
  
  ## --- Step 1: spline derivative term (ζ^T b'_i) ---
  eta_phi <- as.vector(BBD %*% gamma)    # n-vector, ζ^T b'(u_i)
  
  ## --- Step 2: linear predictor and Binomial mean ---
  lin.pred <- BB %*% gamma + X %*% beta
  mu <- plogis(lin.pred)                 # logistic link
  mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)   # numerical stabilization
  
  ## --- Step 3: reparameterization alpha = alpha_eta ---
  theta_alpha <- theta_from_alpha(alpha)     # η vector (length p-1)
  k <- length(theta_alpha)
  r2 <- sum(theta_alpha^2)
  
  ## Normalization factors
  norm1 <- 1 / sqrt(1 + r2)                 # (1+||η||²)^(-1/2)
  norm3 <- norm1 / (1 + r2)                 # (1+||η||²)^(-3/2)
  
  ## --- Step 4: split Z into z1 and z_{-1} ---
  z1    <- Z[, 1]                            # n-vector
  zrest <- Z[, -1, drop = FALSE]             # n × (p-1)
  
  ## Compute z1 + z_{-1}^T η
  zlin <- as.vector(z1 + zrest %*% theta_alpha)
  
  ## --- Step 5: δ_i = derivative of α_η^T z_i wrt η ---
  ## δ_i = norm1 * z_{i,-1}  -  norm3 * (z1_i + z_{i,-1}^T η) * η
  psi_mat <- norm1 * zrest -
    norm3 * (zlin * matrix(theta_alpha,
                           nrow = nrow(Z),
                           ncol = k,
                           byrow = TRUE))   # n × (p-1)
  
  ## --- Step 6: build E-block ---
  ## ξ_i = δ_i * (ζ^T b'(u_i)) = psi_i * eta_phi_i
  E <- rbind(t(psi_mat) %*% diag(eta_phi),   # (p-1) × n  index block
             t(X),                           # q × n
             t(BB))                          # q_n × n
  
  ## --- Step 7: penalty matrix (unchanged) ---
  P <- generate_P_matrix(nrow(E), ncol(BB), 2)
  
  ## --- Step 8: parameter vector (η, β, γ) ---
  theta <- c(theta_alpha, beta, gamma)
  
  #-----------Find log-likelihood function------------
  lik<-sum(Y * log(mu) + (1 - Y) * log(1 - mu)) #-----log likelihood------
  # pen <- 0.5 * lambda * drop(crossprod(theta, P %*% theta))  # (1/2) λ θ^T P θ
  # lik <- sum(Y * log(mu) + (1 - Y) * log(1 - mu)) - pen #-----Penalized log likelihood------
  
  #-----------Find gradient-----------------------
  grad<-E%*%(Y-mu) - lambda*P%*%theta
  
  #-----------Find Hessian matrix------------------
  w  <- pmax(mu * (1 - mu), 1e-8)  # clip to avoid zero weights
  Es <- E * rep(sqrt(w), each = nrow(E))   # scale each column by sqrt(w_i)
  K  <- tcrossprod(Es) + lambda * P        # K = E diag(w) E^T + λP  (SPD target)
  K   <- 0.5 * (K + t(K))        # enforce symmetry
  
  # Try Cholesky on K
  chol_ok <- TRUE
  R <- tryCatch(chol(K), error = function(e) { chol_ok <<- FALSE; NULL })
  
  if (!chol_ok) {
    # Minimal ridge until SPD
    ridge <- 1e-8
    for (iter in 1:8) {
      R <- tryCatch(chol(K + diag(ridge, nrow(K))), error = function(e) NULL)
      if (!is.null(R)) { K <- K + diag(ridge, nrow(K)); break }
      ridge <- ridge * 10
    }
    if (is.null(R)) stop("Failed to make K SPD even after ridge.")
  }
  
  Kinv   <- chol2inv(R)
  H.inv  <- -Kinv       # since H = -K
  
  result<-list(lik,grad,H.inv)
  names(result)[1]<-"lik"
  names(result)[2]<-"grad"
  names(result)[3]<-"H.inv"
  
  return(result)
  
}

FishScore_stephalv_Pois = function(theta,lambda,mfit,maxk=8){
  
  grad.c<-mfit$grad
  H.inv.c<-mfit$H.inv
  for(k in 1:maxk){
    theta.n<-theta-(0.5)^k*H.inv.c%*%grad.c  
    # alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    alpha.n<-alpha_from_theta(theta.n[1:(nAlpha-1)]) # version 2025
    beta.n<-theta.n[nAlpha:(nAlpha+nBeta-1)]
    gamma.n<-theta.n[(nBeta+nAlpha):length(theta.n)]
    new<-findLogLik_pois(Y,X,Z,C,alpha.n,beta.n,gamma.n,lambda)
    lik.diff<-new$lik-mfit$lik
    if(lik.diff>0) break;
  }
  if(lik.diff<=0) theta.n <- theta
  return(theta.n)
  
}


"findLogLik_pois"<-function(Y,X,Z,C,alpha,beta,gamma,lambda){
  
  n = length(Y)
  U<-as.vector(Z%*%alpha)
  dataset<-cbind(Y,X,Z,C,U)
  dataset<-dataset[order(dataset[,ncol(dataset)]),]
  nBeta = ncol(X); nAlpha = ncol(Z)
  Y<-dataset[,1];X<-dataset[,2:(nBeta+1)];Z<-dataset[,(nBeta+2):(nBeta+nAlpha+1)]
  C<-dataset[,(ncol(dataset)-1)];U<-dataset[,ncol(dataset)]
  
  #-----------Find B and derivative of B----------------  
  qn    <- ceiling(n^(1/3))                 ## No. of inner knots
  # qn    <- 4                 ## No. of inner knots
  prob  <- seq(0,1,1/(qn-1))
  psi   <- as.vector(unique(quantile(U,prob)))
  knots <- c(rep(min(U),3),psi,rep(max(U),3))
  B    <- splineDesign(knots = knots,U,ord = 4)  
  BD   <- splineDesign(knots = knots,U,ord = 4, derivs=1)
  
  bs      <- matrix(colSums(B),nrow=1)
  res     <- qr(t(bs))
  bs.Q    <- qr.Q(res,complete = TRUE)
  bs.R    <- qr.R(res,complete = TRUE)
  bs.Q2   <- bs.Q[,seq(ncol(bs.R)+1,ncol(bs.Q))]
  BB   <- B%*%bs.Q2   ## new design matrix compared to original B 
  BBD  <- BD%*%bs.Q2   ## new design matrix compared to original B 
  
  ## eta_phi: derivative of the spline part, plays role of ζ^T b'_i
  eta_phi <- as.vector(BBD %*% gamma)   # n-vector
  
  ## Mean function
  mu <- C * as.vector(exp(BB %*% gamma + X %*% beta))
  
  ## Reparameterization: alpha = alpha_eta, with eta = theta_alpha
  theta_alpha <- theta_from_alpha(alpha)          # eta ∈ R^{p-1}
  k           <- length(theta_alpha)
  r2          <- sum(theta_alpha^2)
  
  ## Normalization factors (1 + ||eta||^2)^(-1/2) and (-3/2)
  norm1 <- 1 / sqrt(1 + r2)                       # (1 + ||eta||^2)^(-1/2)
  norm3 <- norm1 / (1 + r2)                       # (1 + ||eta||^2)^(-3/2)
  
  ## Split Z into first component (z1) and remaining (z_{-1})
  z1    <- Z[, 1]                                 # n-vector
  zrest <- Z[, -1, drop = FALSE]                  # n × (p-1) matrix
  
  ## z1 + z_{-1}^T eta  (n-vector)
  zlin <- as.vector(z1 + zrest %*% theta_alpha)
  
  ## δ_i rows as in the manuscript:
  ## δ_i = (1+r2)^(-1/2) z_{i,-1} - (1+r2)^(-3/2) (z1_i + z_{i,-1}^T eta) * eta
  psi <- norm1 * zrest -
    norm3 * (zlin * matrix(theta_alpha,
                           nrow = nrow(Z),
                           ncol = k,
                           byrow = TRUE))    # n × (p-1)
  
  ## Residuals (Y - µ)
  delta <- Y - mu
  
  ## Index block ξ_i = δ_i * (ζ^T b'_i) implemented as psi^T diag(eta_phi)
  E <- rbind(t(psi) %*% diag(eta_phi),   # (p-1) × n  = index part
             t(X),                       # q × n      = beta part
             t(BB))                      # q_n × n    = gamma (spline) part
  
  ## Penalty matrix (unchanged)
  P <- generate_P_matrix(nrow(E), ncol(BB), 2)
  
  ## Parameter vector in (eta, beta, gamma) order
  theta <- c(theta_alpha, beta, gamma)
  
  #------------------------------------------------
  #
  #  Debug: mu is too large
  #
  #------------------------------------------------
  if(prod(!is.infinite(mu)) == 0){#-----Inf occued
    cat("Debug: mu is infinite!","\n")
    return(-2)
  }
  if(prod(!is.infinite(log(mu))) == 0){#-----Inf occued
    cat("Debug: Log-mu is infinite!","\n")
    mu = mu + 1e-256
  }
  
  #-----------Find log-likelihood function------------
  lik<-sum(Y*log(mu)-mu)#-----log likelihood------
  
  #-----------Find gradient-----------------------
  grad<-E%*%delta - lambda*P%*%theta
  
  
  #-----------Find Hessian matrix------------------
  H<-(E%*%diag(mu)%*%t(E)+(lambda*P))*(-1)
  H <- (H + t(H))/2
  if(prod(!is.infinite(H)) == 0){#-----Inf occued
    cat("Debug: Hessian matrix is infinite!","\n")
    return(-2)
  }
  
  #------------------------------------------------
  #
  #  Debug: Force H to be p.d.
  #
  #------------------------------------------------
  if(!is.positive.definite(-H)){
    cat("Debug: H is NOT positive definte","\n")
    temp<-eigen(-H)$values
    H<-make.positive.definite(-H)*(-1)
  }
  
  #------------------------------------------------
  #
  #  Debug: H is singular
  #
  #------------------------------------------------
  if(rcond(H)<err.rcond){
    cat("Debug: H is singular!","\n")
    H<-H+err.singular
    #return(-1)#-----H is singular
  }
  
  
  
  H.inv<-tryCatch(
    {
      solve(H,diag(ncol(H)))
    },
    error= function(cond){
      message("H is computationally singular")
      return(NULL)
    }
    
  ) 
  
  result<-list(lik,grad,H.inv)
  names(result)[1]<-"lik"
  names(result)[2]<-"grad"
  names(result)[3]<-"H.inv"
  
  return(result)
  
}

# ------------------------------------------------------------------------------
# Alpha reparameterization 
# alpha(theta) = (1, theta_1, ..., theta_{p-1}) / sqrt(1 + ||theta||^2)
# ------------------------------------------------------------------------------

alpha_from_theta <- function(theta) {
  denom <- sqrt(1 + sum(theta^2))
  alpha <- c(1, theta) / denom
  return(alpha)
}

theta_from_alpha <- function(alpha) {
  stopifnot(length(alpha) >= 2)
  alpha1 <- alpha[1]
  if (alpha1 <= 0) warning("alpha1 should be > 0 for identifiability.")
  theta <- alpha[-1] / alpha1
  return(theta)
}

# ------------------------------------------------------------------------------
# Jacobian J = d alpha / d theta   (dimension: p x (p-1))
#   ∂alpha_1/∂theta = -theta / (1 + ||theta||^2)^(3/2)
#   ∂alpha_{j+1}/∂theta_k = delta_{jk}/sqrt(1 + ||theta||^2)
#                           - (theta_j * theta_k)/(1 + ||theta||^2)^(3/2)
# But in compact matrix form:
#   J = rbind(-theta * d^3, d * I_{p-1} - d^3 * theta %*% t(theta))
#   where d = (1 + ||theta||^2)^(-1/2)
# ------------------------------------------------------------------------------

## alpha(theta) = (1, theta) / sqrt(1 + ||theta||^2)
## J = d alpha / d theta  : (p_alpha × p_theta) matrix
jac_alpha_wrt_theta <- function(theta) {
  k <- length(theta)               # p_theta = length(theta)
  s <- sqrt(1 + sum(theta^2))      # s = sqrt(1 + ||theta||^2)
  
  J <- matrix(0, nrow = k + 1, ncol = k)
  
  ## First row: d alpha_1 / d theta = -theta / s^3
  J[1, ] <- -theta / s^3
  
  ## Rows 2..(k+1): d alpha_{j+1} / d theta = e_j / s - theta_j * theta / s^3
  ## This is diag(1/s, k) - (theta %o% theta) / s^3
  J[2:(k + 1), ] <- diag(1 / s, k) - (theta %o% theta) / s^3
  
  return(J)
}


# ------------------------------------------------------------------------------
# Delta-method variance propagation: Var(alpha) = J * Var(theta) * J^T
# V_theta: covariance matrix of unconstrained parameters theta
# ------------------------------------------------------------------------------

## Delta-method: Var(alpha) = J Var(theta) J^T
## theta: length k (theta = alpha_{-1} / alpha_1)
## V_theta: k x k covariance matrix of theta
var_alpha_from_var_theta <- function(theta, V_theta) {
  J <- jac_alpha_wrt_theta(theta)
  V_alpha <- J %*% V_theta %*% t(J)
  return(V_alpha)
}


