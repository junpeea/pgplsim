#-----------------------------------------------
#
#   Single Index Model for Poisson Regresion
#   Recently updated on 11/15/2025
#
rm(list=ls())
# source("sim_lib_251001.R") # 2013 reparametrization
source("sim_lib_251101.R") # 2025 reparametrization
#------------------------------------------------


#------------------------------Main Function-------------------------
n<-1000 #--------Sample size
alpha.true<-c(1,1,1)/sqrt(3)
beta.true<-c(1,-0.5)
eta.true<- alpha.true[-1]/alpha.true[1] #2025 version
phi<- function(t) {
  sin(2*t)
}
nAlpha<-length(alpha.true); nBeta<-length(beta.true)
num_sim = 1000
M<-50
err<-1e-6#---used in step 2
err.rcond<-1e-12#----singular matrix
err.singular<-0.01#---adjust singular matrix

B0<-NULL
nMissing<-nSingular<-nInf<-0

cat("Sample size: ",n,"\n")

sim_results = list()
sim_results$alpha_mean     <- matrix(NA, num_sim, nAlpha)
sim_results$alpha_sd       <- matrix(NA, num_sim, nAlpha)
sim_results$beta_mean      <- matrix(NA, num_sim, nBeta)
sim_results$eta_mean       <- matrix(NA, num_sim, (nAlpha-1))
sim_results$beta_sd        <- matrix(NA, num_sim, nBeta)
sim_results$eta_sd         <- matrix(NA, num_sim, (nAlpha-1))
sim_results$alpha_covprob  <- matrix(NA, num_sim, nAlpha)
sim_results$beta_covprob   <- matrix(NA, num_sim, nBeta)
sim_results$eta_covprob    <- matrix(NA, num_sim, (nAlpha-1))
sim_results$u.h     = list()
sim_results$BB.h    = list()
sim_results$gamma.h = list()
success.ind = rep(FALSE,num_sim)

for(i in 1:num_sim){
  cat("---------------------Iteration ", i,"\n")
  #-----------------Data Generation-----------------
  set.seed(10*i+1); Z1  <- runif(n, 0, 2)     # Predictor Z1
  set.seed(10*i+2); Z2  <- runif(n, 0, 2)     # Predictor Z2
  set.seed(10*i+3); Z3  <- runif(n, 0, 2)     # Predictor Z3
  set.seed(10*i+4); X1  <- rnorm(n, 0, 1)        # Predictor X1
  set.seed(10*i+5); X2  <- rbinom(n, 1, 0.5)     # Predictor X2
  set.seed(10*i+6); C   <-runif(n,0,M)
  
  X<-cbind(X1,X2);Z<-cbind(Z1,Z2,Z3)
  u<- Z%*%alpha.true
  # mu<-exp(X%*%beta.true + phi(u))*C#------Model 1
  mu<-exp(X%*%beta.true + phi(u) - mean(phi(u)))*C#------Mean adjusted model
  Y<-rpois(n,mu)
  
  #---------------Step0: Find initial value of Alpha--------------
  #cat("Step 0","\n")
  if(1){
    result.0<-glm(Y~Z+X, family = poisson(),offset=log(C))$coefficients
    alpha.0<-result.0[2:(nAlpha+1)]
    alpha.0<-alpha.0/sqrt(sum(alpha.0^2))
    if(alpha.0[1]< 0)
      alpha.0<- alpha.0*(-1)
  }else{
    result.0<-glm(Y~-1+Z+X, family = poisson(),offset=log(C))$coefficients
    alpha.0<-result.0[1:nAlpha]
    alpha.0<-alpha.0/sqrt(sum(alpha.0^2))
    if(alpha.0[1]< 0)
      alpha.0<- alpha.0*(-1)
  }
  # alpha.0<- alpha.true
  
  #--------------Step0: Find initial values of Beta and Gamma
  #cat("Step 1","\n")
  u.0<-as.vector(Z%*%alpha.0)
  dataset.0<-cbind(Y,X,Z,C,u.0)
  dataset.1<-dataset.0[order(dataset.0[,ncol(dataset.0)]),]
  Y.1<-dataset.1[,1];X.1<-dataset.1[,2:(nBeta+1)];Z.1<-dataset.1[,(nBeta+2):(nBeta+nAlpha+1)]
  C.1<-dataset.1[,(ncol(dataset.1)-1)];u.1<-dataset.1[,ncol(dataset.1)]
  
  qn    <- ceiling(n^(1/3))                 ## No. of inner knots
  # qn    <- 4
  prob  <- seq(0,1,1/(qn-1))
  psi   <- as.vector(unique(quantile(u.1,prob)))
  knots <- c(rep(min(u.1),3),psi,rep(max(u.1),3))
  B.1   <- splineDesign(knots = knots,u.1,ord = 4) 
  BB.1  <- compute_bmatrix(u.1,B.1,n)
  
  #---------------Plot the graphs--------------------------
  # result.1<-glm.fit(cbind(B.1,X.1), Y.1,family = poisson(),offset=log(C.1))$coefficients
  result.1<-glm.fit(cbind(X.1,BB.1), Y.1,family = poisson(),offset=log(C.1))$coefficients
  
  beta.1<-result.1[1:ncol(X.1)]
  #cat("Beta 0: ", beta.1,"\n")
  gamma.1<-result.1[-c(1:ncol(X.1))]
  
  if(prod(!is.na(gamma.1)) == 0){#----detecting missing values
    nMissing<-nMissing+1
    next
  }
  
  cat("Alpha 0 : ", alpha.0,"\n")
  cat("Beta  0 : ", beta.1,"\n") 
  
  #plot(u.1,u.1,type = "l",col = "red",xlim=c(0,3.5),ylim=c(0,3.5))#---Model 1---
  plot(u.1, phi(u.1)-mean(phi(u.1)),type = "l",col = "brown",lwd = 2)#---Model 1---
  points(u.1,as.vector(BB.1%*%gamma.1),type = "l",col = "blue",lwd = 2,lty = 2)
  
  #--------------Step 1-2: Fisher's scoring approach--------------------------
  alpha.n<-alpha.0; beta.n<-beta.1; gamma.n<-gamma.1
  lik.c<--1;lik.n<-0; lik.diff<-lik.n-lik.c
  lambda.n <- 1
  
  for(iter in 1:M){
    #--------------Step 1: Inner iteration--------------------------
    alpha.c<-alpha.n; beta.c<-beta.n; gamma.c<-gamma.n
    lik.c<-lik.n; lambda.c <- lambda.n
    # theta.c<-c(alpha.c[-1],beta.c,gamma.c) # Version 2013
    theta.c<-c(alpha.c[-1]/alpha.c[1],beta.c,gamma.c) # version 2025
    current<-findLogLik_pois(Y,X,Z,C,alpha.c,beta.c,gamma.c,lambda.c)
    grad.c<-current$grad
    H.inv.c<-current$H.inv
    
    # theta.n<-theta.c-H.inv.c%*%grad.c#-------updated the outcome (step1)
    # theta.n<-FishScore_stephalv_Pois(theta.c,lambda.c,current)#-------updated the outcome (step1)
    
    #--------------Step 2: Outer iteration--------------------------
    lambda.n <- update_lambda(current,alpha.c,beta.c,gamma.c,lambda.c)
    new <-findLogLik_pois(Y,X,Z,C,alpha.c,beta.c,gamma.c,lambda.n)
    grad.n<-new$grad
    H.inv.n<-new$H.inv
    
    # theta.nn<-theta.c-H.inv.n%*%grad.n #-------updated the outcome (step2)
    theta.n<-FishScore_stephalv_Pois(theta.c,lambda.n,new)#-------updated the outcome (step2)
    # alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    # if(sum(is.na(alpha.n))){
    #   alpha.n <- alpha.c
    #   break;
    # }
    # alpha.n<-c(1,theta.n[1:(nAlpha-1)])/sqrt(1+sum(theta.n[1:(nAlpha-1)])^2)  # version 2025
    alpha.n <- alpha_from_theta(theta.n[1:(nAlpha-1)])
    beta.n<-theta.n[nAlpha:(nAlpha+nBeta-1)]
    gamma.n<-theta.n[(nBeta+nAlpha):length(theta.n)]
    
    theta.diff = sqrt(sum((theta.n - theta.c)^2))
    cat("diff: ",theta.diff,"\n")
    
    if(theta.diff<err){
      success.ind[i] <- TRUE
      break
    }
    
  }#while
  alpha.h<-alpha.n; beta.h<-beta.n; gamma.h<-gamma.n
  eta.h<- theta.n[(1:(nAlpha-1))]; lambda.h <- lambda.n
  
  #-----------------Plot the estimated curve------------------------
  u.h<-as.vector(Z%*%alpha.h)
  dataset.h<-cbind(Y,X,Z,C,u.h)
  dataset.h<-dataset.h[order(dataset.h[,ncol(dataset.h)]),]
  Y.h<-dataset.h[,1];X.h<-dataset.h[,2:(nBeta+1)];Z.h<-dataset.h[,(nBeta+2):(nBeta+nAlpha+1)];u.h<-dataset.h[,7]
  C.h<-dataset.h[,(ncol(dataset.h)-1)];u.h<-dataset.h[,ncol(dataset.h)]
  
  psi   <- as.vector(unique(quantile(u.h,prob)))
  knots <- c(rep(min(u.h),3),psi,rep(max(u.h),3))
  B.h   <- splineDesign(knots = knots,u.h,ord = 4)
  BB.h  <- compute_bmatrix(u.h,B.h,n)
  
  cat("Alpha h : ", alpha.h,"\n")
  cat("Beta  h : ", beta.h,"\n")
  
  plot(u.h ,phi(u.h)-mean(phi(u.h)),type = "l",col = "brown",lwd = 2)#---Model 1---
  points(u.h,as.vector(BB.h%*%gamma.h),type = "l",col = "darkgrey",lwd = 2,lty = 5)
  points(u.1,as.vector(BB.1%*%gamma.1),type = "l",col = "blue",lwd = 2,lty = 2)
  
  # ----- Variance estimaion -----
  V_full  <- -H.inv.n
  p_theta <- nAlpha - 1
  V_theta <- V_full[1:p_theta, 1:p_theta]
  
  # ----- Model-based SEs for alpha via delta method on *theta-alpha block only* -----
  V_alpha_model  <- var_alpha_from_var_theta(theta.n[1:p_theta], V_theta)
  se_alpha_model <- sqrt(pmax(diag(V_alpha_model), 0))
  
  # ----- Model-based SEs for beta (standard block) -----
  se_beta_model <- sqrt(pmax(diag(V_full)[nAlpha:(nAlpha + nBeta - 1)], 0))
  se_eta_model <- sqrt(pmax(diag(V_theta), 0))
  
  # store
  sim_results$alpha_mean[i, ] <- alpha.h
  sim_results$beta_mean[i, ]  <- beta.h
  sim_results$eta_mean[i, ]   <- eta.h
  sim_results$alpha_sd[i, ]   <- se_alpha_model
  sim_results$beta_sd[i, ]    <- se_beta_model
  sim_results$eta_sd[i, ]     <- se_eta_model
  
  # coverage using model SEs
  a_lo <- alpha.h + qnorm(0.025) * se_alpha_model
  a_hi <- alpha.h + qnorm(0.975) * se_alpha_model
  sim_results$alpha_covprob[i, ] <- (a_lo < alpha.true) & (a_hi > alpha.true)
  
  b_lo <- beta.h + qnorm(0.025) * se_beta_model
  b_hi <- beta.h + qnorm(0.975) * se_beta_model
  sim_results$beta_covprob[i, ] <- (b_lo < beta.true) & (b_hi > beta.true)
  
  e_lo <- eta.h + qnorm(0.025) * se_eta_model
  e_hi <- eta.h + qnorm(0.975) * se_eta_model
  sim_results$eta_covprob[i, ] <- (e_lo < eta.true) & (e_hi > eta.true)
  
  sim_results$u.h[[i]]     = u.h
  sim_results$BB.h[[i]]    = BB.h
  sim_results$gamma.h[[i]] = gamma.h
  sim_results$lambda.h[[i]] <- lambda.h
  
  cat("Simulation", i, "completed","\n")
}

alpha_mean_pm = colMeans(sim_results$alpha_mean[success.ind,])
beta_mean_pm  = colMeans(sim_results$beta_mean[success.ind,])
eta_mean_pm  = colMeans(sim_results$eta_mean[success.ind,])

# tab1 <-data.frame(
#   Method = c("PM", "PM", "PM", "PM", "PM"),
#   Mean = unlist(c(alpha_mean_pm, beta_mean_pm)),
#   Bias = unlist(c(alpha.true - alpha_mean_pm, beta.true - beta_mean_pm)),
#   SD   = c(apply(sim_results$alpha_mean[success.ind,], 2, sd), apply(sim_results$beta_mean[success.ind,], 2, sd)),
#   mean_of_se = unlist(c(
#     colMeans(sim_results$alpha_sd[success.ind,]),
#     colMeans(sim_results$beta_sd[success.ind,]))),
#   SD_of_se = c(apply(sim_results$alpha_sd[success.ind,], 2, sd),
#                apply(sim_results$beta_sd[success.ind,], 2, sd)), # Standard Deviation of Standard Errors (SD.se)
#   covprob  = c(colMeans(sim_results$alpha_covprob[success.ind,]),
#                colMeans(sim_results$beta_covprob[success.ind,])),
#   MSE = unlist(c(
#     rowMeans(apply(sim_results$alpha_mean[success.ind,], 1, function(v)(v-alpha.true)^2)),
#     rowMeans(apply(sim_results$beta_mean[success.ind,], 1, function(v)(v-beta.true)^2)))),
#   row.names = c("alpha1_pm", "alpha2_pm", "alpha3_pm", "beta1_pm", "beta2_pm")
# )
# tab1

tab2 <-data.frame(
  Mean = unlist(c(alpha_mean_pm, beta_mean_pm,eta_mean_pm)),
  Bias = unlist(c(alpha.true - alpha_mean_pm, beta.true - beta_mean_pm, eta.true - eta_mean_pm)),
  SD   = c(apply(sim_results$alpha_mean[success.ind,], 2, sd), apply(sim_results$beta_mean[success.ind,], 2, sd), apply(sim_results$eta_mean[success.ind,], 2, sd)),
  mean_of_se = unlist(c(
    colMeans(sim_results$alpha_sd[success.ind,]),
    colMeans(sim_results$beta_sd[success.ind,]),
    colMeans(sim_results$eta_sd[success.ind,]))),
  SD_of_se = c(apply(sim_results$alpha_sd[success.ind,], 2, sd),
               apply(sim_results$beta_sd[success.ind,], 2, sd),
               apply(sim_results$eta_sd[success.ind,], 2, sd)), # Standard Deviation of Standard Errors (SD.se)
  covprob  = c(colMeans(sim_results$alpha_covprob[success.ind,]),
               colMeans(sim_results$beta_covprob[success.ind,]),
               colMeans(sim_results$eta_covprob[success.ind,])),
  MSE = unlist(c(
    rowMeans(apply(sim_results$alpha_mean[success.ind,], 1, function(v)(v-alpha.true)^2)),
    rowMeans(apply(sim_results$beta_mean[success.ind,], 1, function(v)(v-beta.true)^2)),
    rowMeans(apply(sim_results$eta_mean[success.ind,], 1, function(v)(v-eta.true)^2)))),
  row.names = c("alpha1_pm", "alpha2_pm", "alpha3_pm", "beta1_pm", "beta2_pm", "eta1_pm", "eta2_pm")
)
tab2

# Mean          Bias          SD  mean_of_se     SD_of_se covprob          MSE
# alpha1_pm  0.5760172  0.0013330936 0.009030783 0.007872653 0.0007067340    0.92 8.251664e-05
# alpha2_pm  0.5776691 -0.0003187887 0.007897659 0.007786242 0.0006624510    0.97 6.185092e-05
# alpha3_pm  0.5781757 -0.0008253953 0.008592053 0.007760101 0.0006031192    0.91 7.376641e-05
# beta1_pm   0.9988957  0.0011042972 0.006909746 0.007198594 0.0008935555    0.95 4.848662e-05
# beta2_pm  -0.4993719 -0.0006280752 0.014259503 0.015045303 0.0006681625    0.98 2.016946e-04
# eta1_pm    1.0032145 -0.0032144840 0.025347783 0.023642490 0.0020485286    0.95 6.464179e-04
# eta2_pm    1.0041298 -0.0041298214 0.027412390 0.023590294 0.0022166436    0.93 7.609802e-04
