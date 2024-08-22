Rcpp::sourceCpp("compute_funcs-20240816.cpp")

#############################################################################
## Generate data for Bernoulli example, create individual and optimal contours 
#############################################################################

expit <- function(z) exp(z)/(1+exp(z))

set.seed(600)
theta <- 0.5
k <- 3
n_i <- c(50,75,100)
theta <- c(-0.5,0.3)
p <- length(theta)
covariates <- lapply(1:k, function(i) cbind(1, rnorm(n_i[i],0,1)) )
Y <- sapply(1:k, function(i) rbinom(n=n_i[i], size=1, prob=expit(covariates[[i]]%*%theta)) )

theta_MLE <- matrix(0,k,p)
J_ni <- list()
for(i in 1:k){
  glm_i <- glm(Y[[i]] ~ 0 + covariates[[i]], family="binomial")
  theta_MLE[i,] <- as.vector(coef(glm_i))
  mu_i <- expit(covariates[[i]]%*%theta_MLE[i,])
  J_ni[[i]] <- Reduce("+", lapply(1:n_i[i], function(j) covariates[[i]][j,] %o% covariates[[i]][j,] * mu_i[j]*(1-mu_i[j])))
}

#range of values for theta
s.J <- Reduce("+",J_ni)
meta_estimate <- drop(solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) ))
theta_vec <- as.matrix(expand.grid(
  seq(meta_estimate[1]-2*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+2*sqrt(diag(solve(s.J)))[1],length=150), 
  seq(meta_estimate[2]-2*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+2*sqrt(diag(solve(s.J)))[2],length=150)) )

# Plotting the individual IM contours
individual_IMs <- lapply(1:k, function(i) {
  valid_logistic_individual_IM(theta=theta_vec, theta_MLE=theta_MLE[i,], Y=Y[[i]], covariates=covariates[[i]], n=n_i[i], p=p, M=50000)
  })

# Plotting the optimal IM contour
combined_glm <- glm(unlist(Y) ~ 0 + do.call(rbind, covariates), family="binomial")
optimal_IM <- valid_logistic_individual_IM(theta=theta_vec, theta_MLE=as.vector(coef(combined_glm)), Y=unlist(Y), 
                                           covariates=do.call(rbind, covariates), n=sum(n_i), p=p, M=50000)

#############################################################################
## Compute Ryan's proposed naive IM
#############################################################################

# function for naive IM
pl.naive <- function(theta, theta_MLE, J_ni, p){
  s.J <- Reduce("+",J_ni)
  theta.n <- solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) )
  apply(theta, 1, function(t) drop(1 - pchisq( t(theta.n -t) %*% s.J %*% (theta.n -t), p)))
}

# Plotting the naive IM
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni, p=p)

#############################################################################
## Compute the proposed middle ground IM with importance weights -- scalar theta, intercept only
#############################################################################

# Middle ground IM
# theta_dag <- as.matrix(expand.grid(
#   c(meta_estimate[1]-sqrt(diag(solve(s.J)))[1], meta_estimate[1], meta_estimate[1]+sqrt(diag(solve(s.J)))[1]), 
#   c(meta_estimate[2]-sqrt(diag(solve(s.J)))[2], meta_estimate[2], meta_estimate[2]+sqrt(diag(solve(s.J)))[2])))
# middle_IM <- logistic_pl_mid(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, 
#                              covariates=covariates, M=50000)

save.image("binomial-settingII-ecdf.RData")

#############################################################################
### Comparison between middle ground IM and "exact" middle ground IM -- scalar theta, intercept only
#############################################################################

middle_IM_exact <- logistic_pl_mid_exact(theta=theta_vec, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, covariates=covariates, M=50000)

save.image("binomial-settingII-ecdf.RData")

#############################################################################
### Empirical CDF of the middle ground contour evaluated at theta_0
#############################################################################

nsim <- 10000
# samples <- <- rep(0,nsim)
exact_samples <- naive_samples <- rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  Y.sim <- sapply(1:k, function(i) rbinom(n=n_i[i], size=1, prob=expit(covariates[[i]]%*%theta)) )
  
  theta_MLE.sim <- matrix(0,k,p)
  J_ni.sim <- list()
  for(i in 1:k){
    glm_i.sim <- glm(Y.sim[[i]] ~ 0 + covariates[[i]], family="binomial")
    theta_MLE.sim[i,] <- as.vector(coef(glm_i.sim))
    mu_i.sim <- expit(covariates[[i]]%*%theta_MLE.sim[i,])
    J_ni.sim[[i]] <- Reduce("+", lapply(1:n_i[i], function(j) covariates[[i]][j,] %o% covariates[[i]][j,] * mu_i.sim[j]*(1-mu_i.sim[j])))
  }
  
  s.J.sim <- Reduce("+",J_ni.sim)
  meta_estimate.sim <- drop(solve(s.J.sim) %*% Reduce("+", lapply(1:length(J_ni.sim), function(i) J_ni.sim[[i]]%*%theta_MLE.sim[i,]) ))
  theta_dag <- meta_estimate.sim + as.matrix(expand.grid(c(-sqrt(diag(solve(s.J.sim)))[1],0,sqrt(diag(solve(s.J.sim)))[1]), 
                                                         c(-sqrt(diag(solve(s.J.sim)))[2],0,sqrt(diag(solve(s.J.sim)))[2])))
  
  # samples[sim] <- logistic_pl_mid(theta=matrix(theta,nrow=1), theta_dag=theta_dag, theta_MLE=theta_MLE.sim, k=k,
  #                                 J_ni=J_ni.sim, n_i=n_i, covariates=covariates, M=50000)
  exact_samples[sim] <- logistic_pl_mid_exact(theta=matrix(theta,nrow=1), theta_MLE=theta_MLE, k=k, 
                                              J_ni=J_ni, n_i=n_i, covariates=covariates, M=50000)
  naive_samples[sim] <- pl.naive(theta=matrix(theta,nrow=1), theta_MLE=theta_MLE.sim, J_ni=J_ni.sim, p=p)
}


save.image("binomial-settingII-ecdf.RData")
