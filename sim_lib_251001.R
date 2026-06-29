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
  
  theta.n<-c(alpha[-1],beta,gamma)
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

FishScore_stephalv_Pois = function(theta,lambda,mfit,maxk=8){
  
  grad.c<-mfit$grad
  H.inv.c<-mfit$H.inv
  for(k in 1:maxk){
    theta.n<-theta-(0.5)^k*H.inv.c%*%grad.c  
    alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    beta.n<-theta.n[nAlpha:(nAlpha+nBeta-1)]
    gamma.n<-theta.n[(nBeta+nAlpha):length(theta.n)]
    new<-findLogLik_pois(Y,X,Z,C,alpha.n,beta.n,gamma.n,lambda)
    lik.diff<-new$lik-mfit$lik
    if(lik.diff>0) break;
  }
  if(lik.diff<=0) theta.n <- theta
  return(theta.n)
  
}

FishScore_stephalv_Bino = function(theta,lambda,mfit,maxk=8){
  
  grad<-mfit$grad
  H.inv<-mfit$H.inv
  for(k in 1:maxk){
    theta.n<-theta-(0.5)^k*H.inv%*%grad  
    alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    if(sum(is.na(alpha.n))){
      alpha.n <- c(0,theta.n[1:(nAlpha-1)])
      break;
    }
    # alpha.n<-c(1,theta.n[1:(nAlpha-1)])/sqrt(1+sum(theta.n[1:(nAlpha-1)]^2))  # version 2025
    beta.n<-theta.n[nAlpha:(nAlpha+nBeta-1)]
    gamma.n<-theta.n[(nBeta+nAlpha):length(theta.n)]
    new<-findLogLik_bino(Y,X,Z,alpha.n,beta.n,gamma.n,lambda)
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
  
  eta<-as.vector(BBD%*%gamma)#------eta------
  mu<-C*as.vector(exp(BB%*%gamma+X%*%beta))#------mu------
  psi<-Z%*%rbind(-alpha[-1]/alpha[1],diag(rep(1,length(alpha[-1]))))#------psi------
  delta<-Y-mu#------delta-------
  E<-rbind(t(psi)%*%diag(eta),t(X),t(BB))#----E----
  P = generate_P_matrix(nrow(E),ncol(BB),2)#----P----
  theta<-c(alpha[-1],beta,gamma)
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
  
  eta<-as.vector(BBD%*%gamma)#------eta------
  lin.pred<- BB%*%gamma+X%*%beta
  mu<-exp(lin.pred)/(1 + exp(lin.pred))#------mu------
  psi<-Z%*%rbind(-alpha[-1]/alpha[1],diag(rep(1,length(alpha[-1]))))#------psi------
  E<-rbind(t(psi)%*%diag(eta),t(X),t(BB))#----E----
  P = generate_P_matrix(nrow(E),ncol(BB),2)#----P----
  theta<-c(alpha[-1],beta,gamma)
  
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
    mu[mu==0] <- 0.000001
  }
  if(prod(!is.infinite(log(1-mu))) == 0){#-----Inf occued
    cat("Debug: Log-mu is infinite!","\n")
    mu[mu==1] <- 0.999999
  }
  
  #-----------Find log-likelihood function------------
  lik<-sum(Y * log(mu) + (1 - Y) * log(1 - mu)) #-----log likelihood------
  
  #-----------Find gradient-----------------------
  grad<-E%*%(Y-mu) - lambda*P%*%theta
  
  
  #-----------Find Hessian matrix------------------
  H<-(E%*%diag(as.vector(mu*(1-mu)))%*%t(E)+(lambda*P))*(-1)
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
