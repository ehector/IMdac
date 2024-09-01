Rcpp::sourceCpp("compute_funcs-20240828.cpp")

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
# individual_IMs <- lapply(1:k, function(i) {
#   valid_logistic_individual_IM(theta=theta_vec, theta_MLE=theta_MLE[i,], Y=Y[[i]], covariates=covariates[[i]], n=n_i[i], p=p, M=50000)
#   })

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
theta_dag <- matrix(meta_estimate, nrow=1)
middle_IM <- logistic_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i,
                             covariates=covariates, M=50000)

save.image("bernoulli-onedag-ecdf.RData")

# dat_optimal_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="oracle", contour=optimal_IM)
dat_difference_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="difference", contour=middle_IM-naive_IM)
dat_naive_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="meta-analysis", contour=naive_IM)
dat_middle_IM <- data.frame(theta_1=theta_vec[,1], theta_2=theta_vec[,2], type="equilibrated", contour=middle_IM)

# dat_IM <- as.data.frame(rbind(dat_optimal_IM, dat_naive_IM, dat_middle_IM))
# dat_IM$type <- factor(dat_IM$type, levels=c("oracle","meta-analysis","equilibrated"), labels=c("oracle","meta-analysis","equilibrated"))
dat_IM <- as.data.frame(rbind(dat_difference_IM, dat_naive_IM, dat_middle_IM))
dat_IM$type <- factor(dat_IM$type, levels=c("meta-analysis","equilibrated","difference"), labels=c("meta-analysis","equilibrated","difference"))
dat_IM$contour[which.max(dat_IM$contour)] <- round(dat_IM$contour[which.max(dat_IM$contour)],2)
dat_IM$contour[dat_IM$type!="difference"][which.min(dat_IM$contour[dat_IM$type!="difference"])] <- 
  round(dat_IM$contour[dat_IM$type!="difference"][which.min(dat_IM$contour[dat_IM$type!="difference"])],1)

plot_all <-
  ggplot(dat_IM[dat_IM$type!="difference",], aes(theta_1, theta_2)) + geom_raster(hjust = 0, vjust = 0, aes(fill = contour)) +
  scale_fill_continuous(type = "viridis", guide=guide_colourbar(title.position="top")) +
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  theme_bw() + theme(axis.title.y=element_text(size=18), axis.title.x=element_text(size=18), legend.title=element_text(size=18),
                     axis.text.y=element_text(size=18), axis.text.x=element_text(size=18),
                     legend.text=element_text(size=18), legend.position="bottom", legend.box = "horizontal", legend.key.width=unit(2,"cm"),
                     strip.background = element_blank(), strip.text.x = element_text(size=18), legend.title.align=0.5) +
  facet_wrap(~type, nrow=1, ncol=2)
plot_naive <-
  ggplot(dat_IM[dat_IM$type=="meta-analysis",], aes(theta_1, theta_2)) + geom_raster(hjust = 0, vjust = 0, aes(fill = contour)) +
  scale_fill_continuous(type = "viridis") +
  xlab("") + ylab(expression(theta[2])) + ggtitle("meta-analysis") +
  theme_bw() + theme(axis.title.y=element_text(size=18), axis.title.x=element_text(size=18), legend.title=element_text(size=18),
                     axis.text.y=element_text(size=18), axis.text.x=element_text(size=18),
                     legend.text=element_text(size=18), plot.title = element_text(size=18, hjust = 0.5))
plot_middle <-
  ggplot(dat_IM[dat_IM$type=="equilibrated",], aes(theta_1, theta_2)) + geom_raster(hjust = 0, vjust = 0, aes(fill = contour)) +
  scale_fill_continuous(type = "viridis") +
  xlab(expression(theta[1])) + ylab("") + ggtitle("equilibrated") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.title.x=element_text(size=18), legend.title=element_text(size=18),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=18),
                     legend.text=element_text(size=18), plot.title=element_text(size=18, hjust = 0.5))
plot_difference <-
  ggplot(dat_IM[dat_IM$type=="difference",], aes(theta_1, theta_2)) + geom_raster(hjust = 0, vjust = 0, aes(fill = contour)) +
  scale_fill_continuous(type = "viridis", option="inferno", guide=guide_colourbar(title.position="top")) +
  xlab("") + ylab("") + ggtitle("difference") + labs(fill = "difference") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.title.x=element_text(size=18), legend.title=element_text(size=18),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=18),
                     legend.key.width=unit(2,"cm"), legend.text=element_text(size=18), legend.position="bottom", legend.box = "horizontal",
                     plot.title = element_text(size=18, hjust = 0.5), legend.title.align=0.5)

legend <- cowplot::get_plot_component(plot_all, 'guide-box-bottom', return_all = TRUE)
legend_diff <- cowplot::get_plot_component(plot_difference, 'guide-box-bottom', return_all = TRUE)
plots_no_legend <- cowplot::plot_grid(
  plot_naive + theme(legend.position="none", plot.margin=margin(5,10,-110,0)),
  plot_middle + theme(legend.position="none", plot.margin=margin(5,10,-113,-10)),
  plot_difference + theme(legend.position="none", plot.margin=margin(5,10,-110,-10)),
  rel_heights=c(1,1,1), rel_widths=c(1.15,1,1), 
  nrow=1, ncol=3)
legends <- cowplot::plot_grid(legend, 
                              legend_diff, 
                              ncol=2, nrow=1, rel_heights=c(1,1), rel_widths=c(2,1))
plots_with_legend <- cowplot::plot_grid(plots_no_legend, legends, ncol=1, nrow=2)

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/logistic-contours.pdf", height=8, width=15)
plots_with_legend
dev.off()

#############################################################################
### Empirical CDF of the middle ground contour evaluated at theta_0
#############################################################################

nsim <- 10000
samples <- naive_samples <- rep(0,nsim)
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
  theta_dag.sim <- as.matrix(expand.grid(
    c(meta_estimate.sim[1]-sqrt(diag(solve(s.J.sim)))[1], meta_estimate.sim[1], meta_estimate.sim[1]+sqrt(diag(solve(s.J.sim)))[1]),
    c(meta_estimate.sim[2]-sqrt(diag(solve(s.J.sim)))[2], meta_estimate.sim[2], meta_estimate.sim[2]+sqrt(diag(solve(s.J.sim)))[2])))
  closest <- which.min(apply(theta_dag.sim, 1, function(x) sum((x-theta)^2) ))
  
  samples[sim] <- logistic_pl_mid_P(theta=matrix(theta,nrow=1), theta_dag=matrix(theta_dag.sim[closest,], nrow=1), theta_MLE=theta_MLE.sim, k=k,
                                  J_ni=J_ni.sim, n_i=n_i, covariates=covariates, M=50000)
  naive_samples[sim] <- pl.naive(theta=matrix(theta,nrow=1), theta_MLE=theta_MLE.sim, J_ni=J_ni.sim, p=p)
}

save.image("bernoulli-settingII-ecdf.RData")

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/ecdf_logistic.pdf", height=5, width=5)
plot(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(samples<=t)),
     type="l", ylab="CDF", xlab=expression(alpha))
lines(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(naive_samples<=t)), col="red")
abline(0,1)
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
dev.off()
