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
n<-400 #--------Sample size
alpha.true<-c(1,1,1)/sqrt(3)
beta.true<-c(1,-0.5)
eta.true<- alpha.true[-1]/alpha.true[1] #2025 version
phi<- function(t) {
  sin(2*t)
}
nAlpha<-length(alpha.true); nBeta<-length(beta.true)
num_sim = 1000
M<-10
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
res = rep(NA, num_sim)

for(i in 1:num_sim){
  
  cat("---------------------Iteration ", i,"\n")
  #-----------------Data Generation-----------------
  set.seed(10*i+1); Z1  <- runif(n, 0, 2)     # Predictor Z1
  set.seed(10*i+2); Z2  <- runif(n, 0, 2)     # Predictor Z2
  set.seed(10*i+3); Z3  <- runif(n, 0, 2)     # Predictor Z3
  set.seed(10*i+4); X1  <- rnorm(n, 0, 1)        # Predictor X1
  set.seed(10*i+5); X2  <- rbinom(n, 1, 0.5)     # Predictor X2
  
  X<-cbind(X1,X2);Z<-cbind(Z1,Z2,Z3)
  u<- Z%*%alpha.true
  lin.pred<- X%*%beta.true + phi(u) - mean(phi(u))#------Mean adjusted model
  mu = exp(lin.pred)/(1 + exp(lin.pred))
  set.seed(10*i+6)
  Y<-rbinom(n, 1, mu)
  res[i] = sum(Y)/length(Y)
  
  #---------------Step0: Find initial value of Alpha--------------
  #cat("Step 0","\n")
  if(1){
    result.0<-glm(Y~Z+X, family = binomial())$coefficients
    alpha.0<-result.0[2:(nAlpha+1)]
    alpha.0<-alpha.0/sqrt(sum(alpha.0^2))
    if(alpha.0[1]< 0)
      alpha.0<- alpha.0*(-1)
  }else{
    result.0<-glm(Y~-1+Z+X, family = binomial())$coefficients
    alpha.0<-result.0[1:nAlpha]
    alpha.0<-alpha.0/sqrt(sum(alpha.0^2))
    if(alpha.0[1]< 0)
      alpha.0<- alpha.0*(-1)
  }
  # alpha.0<- alpha.true
  
  #--------------Step0: Find initial values of Beta and Gamma
  #cat("Step 1","\n")
  u.0<-as.vector(Z%*%alpha.0)
  dataset.0<-cbind(Y,X,Z,u.0)
  dataset.1<-dataset.0[order(dataset.0[,ncol(dataset.0)]),]
  Y.1<-dataset.1[,1];X.1<-dataset.1[,2:(nBeta+1)];Z.1<-dataset.1[,(nBeta+2):(nBeta+nAlpha+1)]
  u.1<-dataset.1[,ncol(dataset.1)]
  
  qn    <- ceiling(n^(1/3))                 ## No. of inner knots
  qn    <- 4
  trim <- 0.01  # remove extreme values
  prob <- seq(trim, 1 - trim, length.out = qn)
  psi   <- as.vector(unique(quantile(u.1,prob)))
  boundary_margin <- 0.05 * (max(u.1) - min(u.1))
  boundary_knots <- c(rep((min(u.1) - boundary_margin),3),psi,rep(max(u.1) + boundary_margin,3))
  B.1   <- splineDesign(knots = boundary_knots,u.1,ord = qn,outer.ok = TRUE) 
  BB.1  <- compute_bmatrix(u.1,B.1,n)
  
  #---------------Plot the graphs--------------------------
  result.1<-glm.fit(cbind(X.1,BB.1), Y.1,family = binomial())$coefficients
  
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
  plot(u.1 ,phi(u.1)-mean(phi(u.1)),type = "l",col = "brown",lwd = 2)
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
    current<-findLogLik_bino(Y,X,Z,alpha.c,beta.c,gamma.c,lambda.c)
    
    grad.c<-current$grad
    H.inv.c<-current$H.inv
    
    # theta.n<-theta.c-H.inv.c%*%grad.c#-------updated the outcome (step1)
    theta.n<-FishScore_stephalv_Bino(theta.c,lambda.c,current)#-------updated the outcome (step1)
    
    #--------------Step 2: Outer iteration--------------------------
    lambda.n <- update_lambda(current,alpha.c,beta.c,gamma.c,lambda.c)
    new <-findLogLik_bino(Y,X,Z,alpha.c,beta.c,gamma.c,lambda.n)
    grad.n<-new$grad
    H.inv.n<-new$H.inv
    
    # theta.nn<-theta.c-H.inv.n%*%grad.n #-------updated the outcome (step2)
    # theta.nn<-FishScore_stephalv_Bino(theta.c,lambda.n,new)#-------updated the outcome (step2)
    # alpha.n<-c(sqrt(1-sum((theta.n[1:(nAlpha-1)])^2)),theta.n[1:(nAlpha-1)])  # version 2013
    # if(sum(is.na(alpha.n))){
    #   alpha.n <- alpha.c
    #   break;
    # }
    alpha.n <- alpha_from_theta(theta.n[1:(nAlpha-1)])  # version 2025
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
  dataset.h<-cbind(Y,X,Z,u.h)
  dataset.h<-dataset.h[order(dataset.h[,ncol(dataset.h)]),]
  Y.h<-dataset.h[,1];X.h<-dataset.h[,2:(nBeta+1)];Z.h<-dataset.h[,(nBeta+2):(nBeta+nAlpha+1)];u.h<-dataset.h[,7]
  u.h<-dataset.h[,ncol(dataset.h)]
  
  psi   <- as.vector(unique(quantile(u.h,prob)))
  boundary_margin <- 0.05 * (max(u.h) - min(u.h))
  boundary_knots <- c(rep((min(u.h) - boundary_margin),3),psi,rep(max(u.h) + boundary_margin,3))
  B.h   <- splineDesign(knots = boundary_knots,u.h,ord = qn,outer.ok = TRUE) 
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

#                 Mean         Bias         SD mean_of_se    SD_of_se   covprob         MSE
# alpha1_pm  0.5712227  0.006127546 0.07565426 0.06667736 0.010118684 0.9080460 0.005755133
# alpha2_pm  0.5728551  0.004495178 0.07674006 0.06654783 0.010097813 0.8996865 0.005903090
# alpha3_pm  0.5727820  0.004568276 0.07660350 0.06655706 0.010394687 0.9049112 0.005882833
# beta1_pm   1.0100173 -0.010017345 0.08825817 0.08713817 0.003547933 0.9540230 0.007881712
# beta2_pm  -0.4965196 -0.003480381 0.10698821 0.10669649 0.003267556 0.9467085 0.011446629
# eta1_pm    1.0335539 -0.033553858 0.27029306 0.21691881 0.087794473 0.9341693 0.074107861
# eta2_pm    1.0330512 -0.033051176 0.27401559 0.21639590 0.090165601 0.9237200 0.076098464