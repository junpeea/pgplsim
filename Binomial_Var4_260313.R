
# source("sim_lib_251101.R") # 2025 reparametrization
source("Updates2026/sim_lib_260305.R") # 2025 reparametrization
source("Updates2026/debug_checklist.R") # 2025 reparametrization

n<-400 #--------Sample size
alpha.true<-c(1,1,1)/sqrt(3)
beta.true<-c(1,-0.5)
eta.true<- alpha.true[-1]/alpha.true[1] #2025 version
phi<- function(t) {
  sin(2*t)
}
nAlpha<-length(alpha.true)
nBeta<-length(beta.true)
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
dbg_list = list()
for(i in 1:num_sim){
  
  #  cat("---------------------Iteration ", i,"\n")
  
  #-----------------Data Generation-----------------
  set.seed(10*i+1) 
  Z1  <- runif(n, 0, 2)     # Predictor Z1
  
  set.seed(10*i+2)
  Z2  <- runif(n, 0, 2)     # Predictor Z2
  
  set.seed(10*i+3) 
  Z3  <- runif(n, 0, 2)     # Predictor Z3
  
  set.seed(10*i+4) 
  X1  <- rnorm(n, 0, 1)     # Predictor X1
  
  set.seed(10*i+5); 
  X2  <- rbinom(n, 1, 0.5)  # Predictor X2
  
  X<-cbind(X1,X2)
  Z<-cbind(Z1,Z2,Z3)
  
  u<- Z%*%alpha.true
  lin.pred<- X%*%beta.true + phi(u) - mean(phi(u)) #------Mean adjusted model
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
  Y.1<-dataset.1[,1]
  X.1<-dataset.1[,2:(nBeta+1)]
  Z.1<-dataset.1[,(nBeta+2):(nBeta+nAlpha+1)]
  u.1<-dataset.1[,ncol(dataset.1)]
  
  knot_spec <- list(
    ord = 4,                  # splineDesign ord (e.g., 4 for cubic)
    qn  = ceiling(n^(1/4)),         ## No. of inner knots
    pen_order = 2              # 2nd-difference penalty
  )
  trim <- 0.01  # remove extreme values
  prob <- seq(trim, 1 - trim, length.out = knot_spec$qn)
  psi   <- as.vector(unique(quantile(u.1,prob)))
  boundary_margin <- 0.05 * (max(u.1) - min(u.1))
  knot_spec$knots <- c(rep((min(u.1) - boundary_margin),(knot_spec$ord-1)),psi,rep(max(u.1) + boundary_margin,(knot_spec$ord-1)))
  B.1   <- splineDesign(knots = knot_spec$knots,u.1,ord = knot_spec$ord,outer.ok = TRUE)
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
  
  #  cat("Alpha 0 : ", alpha.0,"\n")
  #  cat("Beta  0 : ", beta.1,"\n") 
  
  #  plot(u.1,u.1,type = "l",col = "red",xlim=c(0,3.5),ylim=c(0,3.5))#---Model 1---
  #  plot(u.1 ,phi(u.1)-mean(phi(u.1)),type = "l",col = "brown",lwd = 2)
  #  points(u.1,as.vector(BB.1%*%gamma.1),type = "l",col = "blue",lwd = 2,lty = 2)
  
  
  #--------------Step 1-2: Fisher's scoring approach--------------------------
  alpha.n<-alpha.0 
  beta.n<-beta.1 
  gamma.n<-gamma.1
  
  lik.c<--1
  lik.n<-0
  lik.diff<-lik.n-lik.c
  lambda.n <- 1
  
  for(iter in 1:M){
    #--------------Step 1: Inner iteration--------------------------
    alpha.c<-alpha.n; beta.c<-beta.n; gamma.c<-gamma.n
    lik.c<-lik.n; lambda.c <- lambda.n
    # theta.c<-c(alpha.c[-1],beta.c,gamma.c) # Version 2013
    theta.c<-c(alpha.c[-1]/alpha.c[1],beta.c,gamma.c) # version 2025
    current<-findLogLik_bino(Y,X,Z,alpha.c,beta.c,gamma.c,lambda.c,knot_spec)
    
    grad.c<-current$grad
    H.inv.c<-current$H.inv
    
    # theta.n<-theta.c-H.inv.c%*%grad.c#-------updated the outcome (step1)
    theta.n<-FishScore_stephalv_Bino(theta.c,lambda.c,current,knot_spec,Y,X,Z)#-------updated the outcome (step1)
    
    #--------------Step 2: Outer iteration--------------------------
    lambda.n <- update_lambda(current,alpha.c,beta.c,gamma.c,lambda.c)
    new <-findLogLik_bino(Y,X,Z,alpha.c,beta.c,gamma.c,lambda.n,knot_spec)
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
    #   cat("diff: ",theta.diff,"\n")
    
    if(theta.diff<err){
      success.ind[i] <- TRUE
      break
    }
    
  }#while
  alpha.h<-alpha.n 
  beta.h<-beta.n 
  gamma.h<-gamma.n
  eta.h<- theta.n[(1:(nAlpha-1))] 
  lambda.h <- lambda.n
  
  #-----------------Plot the estimated curve------------------------
  u.h<-as.vector(Z%*%alpha.h)
  dataset.h<-cbind(Y,X,Z,u.h)
  dataset.h<-dataset.h[order(dataset.h[,ncol(dataset.h)]),]
  Y.h<-dataset.h[,1]
  X.h<-dataset.h[,2:(nBeta+1)]
  Z.h<-dataset.h[,(nBeta+2):(nBeta+nAlpha+1)]
  u.h<-dataset.h[,7]
  u.h<-dataset.h[,ncol(dataset.h)]
  
  B.h   <- splineDesign(knots = knot_spec$knots,u.h,ord = knot_spec$ord,outer.ok = TRUE) 
  BB.h  <- compute_bmatrix(u.h,B.h,n)
  
  cat("Alpha h : ", alpha.h,"\n")
  cat("Beta  h : ", beta.h,"\n")
  
  plot(u.h ,phi(u.h)-mean(phi(u.h)),type = "l",col = "brown",lwd = 2)#---Model 1---
  points(u.h,as.vector(BB.h%*%gamma.h),type = "l",col = "darkgrey",lwd = 2,lty = 5)
  points(u.1,as.vector(BB.1%*%gamma.1),type = "l",col = "blue",lwd = 2,lty = 2)
  
  # ============================================================
  # Variance estimation (4 methods): model / robust / leverage / finite-sample
  #   - for logistic penalized single-index spline fit
  #   - matches Section 4 of your manuscript
  # ============================================================
  
  # ----------------- run var estimators at final fit -----------------
  var4 <- var4_bino(
    Y = Y, X = X, Z = Z,
    alpha = alpha.h, beta = beta.h, gamma = gamma.h,
    lambda = lambda.h, knot_spec=knot_spec
  )
  
  # unpack
  V_full_model  <- var4$V_model
  V_full_robust <- var4$V_robust
  V_full_lev    <- var4$V_levadj
  V_full_fs     <- var4$V_finite
  
  # ---- blocks for eta/beta/alpha ----
  p_theta <- nAlpha - 1
  idx_eta   <- 1:p_theta
  idx_beta  <- nAlpha:(nAlpha + nBeta - 1)
  
  V_theta_model  <- V_full_model [idx_eta,  idx_eta,  drop = FALSE]
  V_theta_robust <- V_full_robust[idx_eta,  idx_eta,  drop = FALSE]
  V_theta_lev    <- V_full_lev   [idx_eta,  idx_eta,  drop = FALSE]
  V_theta_fs     <- V_full_fs    [idx_eta,  idx_eta,  drop = FALSE]
  
  # alpha via delta method (you already have var_alpha_from_var_theta())
  V_alpha_model  <- var_alpha_from_var_theta(theta.n[idx_eta], V_theta_model)
  V_alpha_robust <- var_alpha_from_var_theta(theta.n[idx_eta], V_theta_robust)
  V_alpha_lev    <- var_alpha_from_var_theta(theta.n[idx_eta], V_theta_lev)
  V_alpha_fs     <- var_alpha_from_var_theta(theta.n[idx_eta], V_theta_fs)
  
  se_alpha_model  <- sqrt(pmax(diag(V_alpha_model),  0))
  se_alpha_robust <- sqrt(pmax(diag(V_alpha_robust), 0))
  se_alpha_lev    <- sqrt(pmax(diag(V_alpha_lev),    0))
  se_alpha_fs     <- sqrt(pmax(diag(V_alpha_fs),     0))
  
  # beta SEs (block diagonal from full V)
  se_beta_model  <- sqrt(pmax(diag(V_full_model )[idx_beta],  0))
  se_beta_robust <- sqrt(pmax(diag(V_full_robust)[idx_beta],  0))
  se_beta_lev    <- sqrt(pmax(diag(V_full_lev   )[idx_beta],  0))
  se_beta_fs     <- sqrt(pmax(diag(V_full_fs    )[idx_beta],  0))
  
  # eta SEs (direct)
  se_eta_model  <- sqrt(pmax(diag(V_theta_model),  0))
  se_eta_robust <- sqrt(pmax(diag(V_theta_robust), 0))
  se_eta_lev    <- sqrt(pmax(diag(V_theta_lev),    0))
  se_eta_fs     <- sqrt(pmax(diag(V_theta_fs),     0))
  
  # ----------------- store (non-breaking: keep your existing fields) -----------------
  # Keep your original model-based storage:
  sim_results$alpha_mean[i, ] <- alpha.h
  sim_results$beta_mean[i, ]  <- beta.h
  sim_results$eta_mean[i, ]   <- eta.h
  sim_results$alpha_sd[i, ] <- se_alpha_model
  sim_results$beta_sd[i, ]  <- se_beta_model
  sim_results$eta_sd[i, ]   <- se_eta_model
  
  # Add new storage for other 3 methods (create on the fly)
  if (is.null(sim_results$alpha_sd_robust)) sim_results$alpha_sd_robust <- sim_results$alpha_sd
  if (is.null(sim_results$beta_sd_robust))  sim_results$beta_sd_robust  <- sim_results$beta_sd
  if (is.null(sim_results$eta_sd_robust))   sim_results$eta_sd_robust   <- sim_results$eta_sd
  
  if (is.null(sim_results$alpha_sd_levadj)) sim_results$alpha_sd_levadj <- sim_results$alpha_sd
  if (is.null(sim_results$beta_sd_levadj))  sim_results$beta_sd_levadj  <- sim_results$beta_sd
  if (is.null(sim_results$eta_sd_levadj))   sim_results$eta_sd_levadj   <- sim_results$eta_sd
  
  if (is.null(sim_results$alpha_sd_finite)) sim_results$alpha_sd_finite <- sim_results$alpha_sd
  if (is.null(sim_results$beta_sd_finite))  sim_results$beta_sd_finite  <- sim_results$beta_sd
  if (is.null(sim_results$eta_sd_finite))   sim_results$eta_sd_finite   <- sim_results$eta_sd
  
  sim_results$alpha_sd_robust[i, ] <- se_alpha_robust
  sim_results$beta_sd_robust[i, ]  <- se_beta_robust
  sim_results$eta_sd_robust[i, ]   <- se_eta_robust
  
  sim_results$alpha_sd_levadj[i, ] <- se_alpha_lev
  sim_results$beta_sd_levadj[i, ]  <- se_beta_lev
  sim_results$eta_sd_levadj[i, ]   <- se_eta_lev
  
  sim_results$alpha_sd_finite[i, ] <- se_alpha_fs
  sim_results$beta_sd_finite[i, ]  <- se_beta_fs
  sim_results$eta_sd_finite[i, ]   <- se_eta_fs
  
  # ----------------- coverage for all 4 methods (optional but useful) -----------------
  # (1) model-based (your existing)
  a_lo <- alpha.h + qnorm(0.025) * se_alpha_model
  a_hi <- alpha.h + qnorm(0.975) * se_alpha_model
  sim_results$alpha_covprob[i, ] <- (a_lo < alpha.true) & (a_hi > alpha.true)
  
  b_lo <- beta.h + qnorm(0.025) * se_beta_model
  b_hi <- beta.h + qnorm(0.975) * se_beta_model
  sim_results$beta_covprob[i, ] <- (b_lo < beta.true) & (b_hi > beta.true)
  
  e_lo <- eta.h + qnorm(0.025) * se_eta_model
  e_hi <- eta.h + qnorm(0.975) * se_eta_model
  sim_results$eta_covprob[i, ] <- (e_lo < eta.true) & (e_hi > eta.true)
  
  # create additional covprob slots if missing
  if (is.null(sim_results$alpha_covprob_robust)) sim_results$alpha_covprob_robust <- sim_results$alpha_covprob
  if (is.null(sim_results$beta_covprob_robust)) sim_results$beta_covprob_robust <- sim_results$beta_covprob
  if (is.null(sim_results$eta_covprob_robust)) sim_results$eta_covprob_robust <- sim_results$eta_covprob
  if (is.null(sim_results$alpha_covprob_levadj)) sim_results$alpha_covprob_levadj <- sim_results$alpha_covprob
  if (is.null(sim_results$beta_covprob_levadj)) sim_results$beta_covprob_levadj <- sim_results$beta_covprob
  if (is.null(sim_results$eta_covprob_levadj)) sim_results$eta_covprob_levadj <- sim_results$eta_covprob
  if (is.null(sim_results$alpha_covprob_finite)) sim_results$alpha_covprob_finite <- sim_results$alpha_covprob
  if (is.null(sim_results$beta_covprob_finite)) sim_results$beta_covprob_finite <- sim_results$beta_covprob
  if (is.null(sim_results$eta_covprob_finite)) sim_results$eta_covprob_finite <- sim_results$eta_covprob
  
  # robust
  a_lo <- alpha.h + qnorm(0.025) * se_alpha_robust
  a_hi <- alpha.h + qnorm(0.975) * se_alpha_robust
  sim_results$alpha_covprob_robust[i, ] <- (a_lo < alpha.true) & (a_hi > alpha.true)
  
  b_lo <- beta.h + qnorm(0.025) * se_beta_robust
  b_hi <- beta.h + qnorm(0.975) * se_beta_robust
  sim_results$beta_covprob_robust[i, ] <- (b_lo < beta.true) & (b_hi > beta.true)
  
  e_lo <- eta.h + qnorm(0.025) * se_eta_robust
  e_hi <- eta.h + qnorm(0.975) * se_eta_robust
  sim_results$eta_covprob_robust[i, ] <- (e_lo < eta.true) & (e_hi > eta.true)
  
  # leverage-adjusted
  a_lo <- alpha.h + qnorm(0.025) * se_alpha_lev
  a_hi <- alpha.h + qnorm(0.975) * se_alpha_lev
  sim_results$alpha_covprob_levadj[i, ] <- (a_lo < alpha.true) & (a_hi > alpha.true)
  
  b_lo <- beta.h + qnorm(0.025) * se_beta_lev
  b_hi <- beta.h + qnorm(0.975) * se_beta_lev
  sim_results$beta_covprob_levadj[i, ] <- (b_lo < beta.true) & (b_hi > beta.true)
  
  e_lo <- eta.h + qnorm(0.025) * se_eta_lev
  e_hi <- eta.h + qnorm(0.975) * se_eta_lev
  sim_results$eta_covprob_levadj[i, ] <- (e_lo < eta.true) & (e_hi > eta.true)
  
  # finite-sample corrected
  a_lo <- alpha.h + qnorm(0.025) * se_alpha_fs
  a_hi <- alpha.h + qnorm(0.975) * se_alpha_fs
  sim_results$alpha_covprob_finite[i, ] <- (a_lo < alpha.true) & (a_hi > alpha.true)
  
  b_lo <- beta.h + qnorm(0.025) * se_beta_fs
  b_hi <- beta.h + qnorm(0.975) * se_beta_fs
  sim_results$beta_covprob_finite[i, ] <- (b_lo < beta.true) & (b_hi > beta.true)
  
  e_lo <- eta.h + qnorm(0.025) * se_eta_fs
  e_hi <- eta.h + qnorm(0.975) * se_eta_fs
  sim_results$eta_covprob_finite[i, ] <- (e_lo < eta.true) & (e_hi > eta.true)
  
  # (optional) keep diagnostics
  if (is.null(sim_results$var_extras)) sim_results$var_extras <- vector("list", num_sim)
  sim_results$var_extras[[i]] <- var4$extras
  
  sim_results$u.h[[i]]     = u.h
  sim_results$BB.h[[i]]    = BB.h
  sim_results$gamma.h[[i]] = gamma.h
  sim_results$lambda.h[[i]] <- lambda.h
  
  cat("Simulation", i, "completed","\n")
  
  dbg_list[[i]] <- debug_checklist_bino(Y, X, Z,
                                        alpha.h, beta.h, gamma.h, lambda.h,
                                        knot_spec,
                                        clamp_h = 0.99)
  
}  

alpha_mean_pm = colMeans(sim_results$alpha_mean)
beta_mean_pm  = colMeans(sim_results$beta_mean)
eta_mean_pm  = colMeans(sim_results$eta_mean)

tab2 <-data.frame(
  Mean = unlist(c(alpha_mean_pm, beta_mean_pm,eta_mean_pm)),
  Bias = unlist(c(alpha.true - alpha_mean_pm, beta.true - beta_mean_pm, eta.true - eta_mean_pm)),
  SD   = c(apply(sim_results$alpha_mean, 2, sd), apply(sim_results$beta_mean, 2, sd), apply(sim_results$eta_mean, 2, sd)),
  SE_mod = unlist(c(
    colMeans(sim_results$alpha_sd),
    colMeans(sim_results$beta_sd),
    colMeans(sim_results$eta_sd))),
  SE_rob = unlist(c(
    colMeans(sim_results$alpha_sd_robust),
    colMeans(sim_results$beta_sd_robust),
    colMeans(sim_results$eta_sd_robust))),
  SE_lev = unlist(c(
    colMeans(sim_results$alpha_sd_levadj),
    colMeans(sim_results$beta_sd_levadj),
    colMeans(sim_results$eta_sd_levadj))),
  SE_fin = unlist(c(
    colMeans(sim_results$alpha_sd_finite),
    colMeans(sim_results$beta_sd_finite),
    colMeans(sim_results$eta_sd_finite))),
  covp_mod  = c(colMeans(sim_results$alpha_covprob),
                colMeans(sim_results$beta_covprob),
                colMeans(sim_results$eta_covprob)),
  covp_rob  = c(colMeans(sim_results$alpha_covprob_robust),
                colMeans(sim_results$beta_covprob_robust),
                colMeans(sim_results$eta_covprob_robust)),
  covp_lev  = c(colMeans(sim_results$alpha_covprob_levadj),
                colMeans(sim_results$beta_covprob_levadj),
                colMeans(sim_results$eta_covprob_levadj)),
  covp_fin  = c(colMeans(sim_results$alpha_covprob_finite),
                colMeans(sim_results$beta_covprob_finite),
                colMeans(sim_results$eta_covprob_finite)),
  MSE = unlist(c(
    rowMeans(apply(sim_results$alpha_mean, 1, function(v)(v-alpha.true)^2)),
    rowMeans(apply(sim_results$beta_mean, 1, function(v)(v-beta.true)^2)),
    rowMeans(apply(sim_results$eta_mean, 1, function(v)(v-eta.true)^2)))),
  row.names = c("alpha1_pm", "alpha2_pm", "alpha3_pm", "beta1_pm", "beta2_pm", "eta1_pm", "eta2_pm")
)
dbg_df <- do.call(rbind, lapply(dbg_list, as.data.frame))
summary(dbg_df)

tab2
save.image(file=paste0("pm binomial ",n,"n var4_260421.RData"))



