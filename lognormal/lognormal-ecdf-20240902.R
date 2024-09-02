Rcpp::sourceCpp("compute_funcs-20240828.cpp")
# library(ggplot2)

#############################################################################
## Generate data for logNormal example, create individual and optimal contours
#############################################################################

set.seed(600)
n_i <- c(12,14,10,11,13)
k <- length(n_i)
mu <- 1
gamma <- 1.5
theta <- c(mu, gamma)
p <- length(theta)
Y <- lapply(1:k, function(i) rlnorm(n=n_i[i], meanlog=mu, sdlog=sqrt(gamma)) )

theta_MLE <- t(sapply(1:k, function(i) lognormal_MLE(Y[[i]]) ))
J_ni <- lapply(1:k, function(i) lognormal_FI(Y[[i]], n_i[i], theta_MLE[i,]))

#range of values for theta
s.J <- Reduce("+",J_ni)
meta_estimate <- drop(solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) ))
theta_vec <- as.matrix(expand.grid(
  seq(meta_estimate[1]-5*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+5*sqrt(diag(solve(s.J)))[1],length=100), 
  seq(meta_estimate[2]-5*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+5*sqrt(diag(solve(s.J)))[2],length=100)) )

# Individual IM contours
individual_IMs <- lapply(1:k, function(i) valid_lognormal_individual_IM(theta=theta_vec, Y=Y[[i]], n=n_i[i], M=50000) )

# Optimal IM contour
optimal_IM <- valid_lognormal_individual_IM(theta=theta_vec, Y=unlist(Y), n=sum(n_i), M=50000)

#############################################################################
## Compute Ryan's proposed naive IM
#############################################################################

# function for naive IM
pl.naive <- function(theta, theta_MLE, J_ni, p){
  s.J <- Reduce("+",J_ni)
  theta.n <- solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) )
  apply(theta, 1, function(t) drop(1 - pchisq( t(theta.n -t) %*% s.J %*% (theta.n -t), p)))
}

# Naive IM
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni, p=p)

#############################################################################
### Comparison between middle ground IM and "exact" middle ground IM
#############################################################################

# Middle ground IM
theta_dag <- matrix(meta_estimate, nrow=1)
middle_IM <- lognormal_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i=n_i, M=50000)

dat_optimal_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="oracle", contour=optimal_IM)
dat_naive_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="meta-analysis", contour=naive_IM)
dat_middle_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="equilibrated", contour=middle_IM)

dat_IM <- as.data.frame(rbind(dat_optimal_IM, dat_naive_IM, dat_middle_IM))
dat_IM$type <- factor(dat_IM$type, levels=c("oracle","meta-analysis","equilibrated"), labels=c("oracle","meta-analysis","equilibrated"))
dat_IM$contour[which.max(dat_IM$contour)] <- round(dat_IM$contour[which.max(dat_IM$contour)],2)

plot_IMs <-
  ggplot(dat_IM, aes(theta_1, theta_2)) + geom_raster(hjust = 0, vjust = 0, aes(fill = contour)) +
  scale_fill_continuous(type = "viridis") +
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  theme_bw() + theme(axis.title.y=element_text(size=18), axis.title.x=element_text(size=18), legend.title=element_text(size=18),
                     axis.text.y=element_text(size=18), axis.text.x=element_text(size=18),
                     legend.text=element_text(size=18),
                     strip.background = element_blank(), strip.text.x = element_text(size=18)) +
  facet_wrap(~type, nrow=1, ncol=3)

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/lognormal-contours.pdf", height=5, width=15)
plot_IMs
dev.off()

#############################################################################
### Empirical CDF of the middle ground contour evaluated at theta_0
#############################################################################

nsim <- 10000
samples <- naive_samples <- rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  Y.sim <- lapply(1:k, function(i) rlnorm(n=n_i[i], meanlog=mu, sdlog=sqrt(gamma)) )
  theta_MLE.sim <- t(sapply(1:k, function(i) lognormal_MLE(Y.sim[[i]]) ))
  J_ni.sim <- lapply(1:k, function(i) lognormal_FI(Y.sim[[i]], n_i[i], theta_MLE.sim[i,]))
  
  s.J.sim <- Reduce("+",J_ni.sim)
  meta_estimate.sim <- drop(solve(s.J.sim) %*% Reduce("+", lapply(1:length(J_ni.sim), function(i) J_ni.sim[[i]]%*%theta_MLE.sim[i,]) ))
  
  theta_dag.sim <- matrix(meta_estimate.sim, nrow=1)
  
  samples[sim] <- lognormal_pl_mid_P(theta=matrix(theta,nrow=1), theta_dag=theta_dag.sim, theta_MLE=theta_MLE.sim, k=k, 
                                     J_ni=J_ni.sim, n_i, M=50000)
  naive_samples[sim] <- pl.naive(theta=matrix(theta,nrow=1), theta_MLE=theta_MLE.sim, J_ni=J_ni.sim, p=p)
}

save.image("lognormal-ecdf.RData")

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/ecdf_lognormal.pdf", height=5, width=5)
plot(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(samples<=t)),
     type="l", ylab="CDF", xlab=expression(alpha))
lines(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(naive_samples<=t)), col="red")
abline(0,1)
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
dev.off()
