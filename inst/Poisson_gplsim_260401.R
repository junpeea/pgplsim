
# source("sim_lib_251101.R") # 2025 reparametrization
source("Updates2026/sim_lib_260220.R") # 2025 reparametrization
source("Updates2026/debug_checklist.R") # 2025 reparametrization
library(gplsim)

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
cat("Sample size: ",n,"\n")

sim_results = list()
sim_results$alpha_mean     <- matrix(NA, num_sim, nAlpha)
sim_results$beta_mean      <- matrix(NA, num_sim, nBeta)
success.ind = rep(FALSE,num_sim)
res = rep(NA, num_sim)

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
  
  set.seed(10*i+6); C   <-runif(n,0,M)    # Offset
  
  X<-cbind(X1,X2)
  Z<-cbind(Z1,Z2,Z3)
  
  u<- Z%*%alpha.true
  mu<-exp(X%*%beta.true + phi(u) - mean(phi(u)))*C#------Mean adjusted model
  Y<-rpois(n,mu)
  
  gpl_fit <- gplsim(
    Y = Y, 
    X = Z, 
    Z = X, 
    family = poisson(link = "log"),  # Update to logit for binary outcomes 
  )
  names(gpl_fit)
  sim_results$alpha_mean[i, ] <- gpl_fit$theta
  sim_results$beta_mean[i, ]  <- gpl_fit$gamma
  
  cat("Simulation", i, "completed",paste0(Sys.time()),"\n")
}  

alpha_mean_Yu = colMeans(sim_results$alpha_mean)
beta_mean_Yu  = colMeans(sim_results$beta_mean)


tab2 <-data.frame(
  Mean = unlist(c(alpha_mean_Yu, beta_mean_Yu)),
  Bias = unlist(c(alpha.true - alpha_mean_Yu, beta.true - beta_mean_Yu)),
  SD   = c(apply(sim_results$alpha_mean, 2, sd), apply(sim_results$beta_mean, 2, sd)),
  MSE = unlist(c(
    rowMeans(apply(sim_results$alpha_mean, 1, function(v)(v-alpha.true)^2)),
    rowMeans(apply(sim_results$beta_mean, 1, function(v)(v-beta.true)^2)))),
  row.names = c("alpha1_gpl", "alpha2_gpl", "alpha3_gpl", "beta1_gpl", "beta2_gpl")
)
tab2


