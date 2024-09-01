Rcpp::sourceCpp("compute_funcs-20240828.cpp")

#############################################################################
## Generate data for Bernoulli example, create individual and optimal contours -- scalar theta, intercept only
#############################################################################

expit <- function(z) exp(z)/(1+exp(z))

sim_vec <- setdiff(1:504, c(32,51,224,319))
sim <- sim_vec[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

set.seed(sim)
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
  seq(meta_estimate[1]-2*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+2*sqrt(diag(solve(s.J)))[1],length=400), 
  seq(meta_estimate[2]-2*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+2*sqrt(diag(solve(s.J)))[2],length=400)) )

theta_vec <- rbind(theta, theta_vec)

# Individual IM contours
# individual_IMs <- lapply(1:k, function(i) {
#   valid_logistic_individual_IM(theta=theta_vec, theta_MLE=theta_MLE[i,], Y=Y[[i]], covariates=covariates[[i]], n=n_i[i], p=p, M=50000) 
#   })
# 
# Optimal IM contour
# combined_glm <- glm(unlist(Y) ~ 0 + do.call(rbind, covariates), family="binomial")
# optimal_IM <- valid_logistic_individual_IM(theta=theta_vec, theta_MLE=as.vector(coef(combined_glm)), Y=unlist(Y), 
#                                            covariates=do.call(rbind, covariates), n=sum(n_i), p=p, M=50000)

#############################################################################
## Compute Ryan's proposed naive IM
#############################################################################

# function for naive IM
pl.naive <- function(theta, theta_MLE, J_ni){
  s.J <- Reduce("+",J_ni)
  theta.n <- solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) )
  apply(theta, 1, function(t) drop(1 - pchisq( t(theta.n -t) %*% s.J %*% (theta.n -t), 1)))
}

# Plotting the naive IM
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni)

#############################################################################
## Compute the proposed middle ground IM with importance weights -- vector theta
#############################################################################

# Plotting middle ground IM
theta_dag <- as.matrix(expand.grid(
  c(meta_estimate[1]-sqrt(diag(solve(s.J)))[1], meta_estimate[1], meta_estimate[1]+sqrt(diag(solve(s.J)))[1]),
  c(meta_estimate[2]-sqrt(diag(solve(s.J)))[2], meta_estimate[2], meta_estimate[2]+sqrt(diag(solve(s.J)))[2])))
middle_IM <- logistic_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i,
                               covariates=covariates, M=50000)

save.image(paste0("bernoulli-settingII-sim",sim,".RData"))
