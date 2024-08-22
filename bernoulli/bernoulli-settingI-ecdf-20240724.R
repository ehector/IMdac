#############################################################################
## Generate data for Bernoulli example, create individual and optimal contours -- scalar theta, intercept only
#############################################################################

expit <- function(z) exp(z)/(1+exp(z))

set.seed(600)
theta <- 0.5
k <- 3
n_i <- c(10,15,20)
theta <- -0.5
Y <- sapply(1:k, function(i) rbinom(n=1, size=n_i[i], prob=expit(theta)) )

theta_MLE <- log(Y/n_i/(1-Y/n_i))
J_ni <- n_i * exp(theta_MLE) / (1+exp(theta_MLE))^2

#range of values for theta
theta_vec <- seq(-3,1,0.01)

logistic_individual_IM <- function(theta, Y, n){
  dbinom(Y, size=n, prob=expit(theta)) / dbinom(Y, size=n, prob=Y/n)
}
valid_logistic_individual_IM <- function(theta, Y, n, M){
  observed_contour <- logistic_individual_IM(theta=theta, Y=Y, n=n)
  Ysim <- sapply(theta, function(t) sapply(1:M, function(sim) rbinom(1, size=n, prob=expit(t)) ))
  pl <- rep(0, length(theta))
  for(t in 1:length(theta)){
    res <- sapply(1:M, function(sim) {
      if(Ysim[sim,t]==0) theta_MLE_sim <- -Inf
      else if (Ysim[sim,t]==n) theta_MLE_sim <- Inf
      else theta_MLE_sim <- log(Ysim[sim,t]/n/(1-Ysim[sim,t]/n))
      
      as.numeric( logistic_individual_IM(theta=theta[t], Y=Ysim[sim,t], n=n) <= observed_contour[t])
    })
    pl[t] <-   mean(res)
  }
  return(pl)
}

# Plotting the individual IM contours
individual_IMs <- lapply(1:k, function(i) valid_logistic_individual_IM(theta=theta_vec, Y=Y[i], n=n_i[i], M=50000) )
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")

# Plotting the optimal IM contour
optimal_IM <- valid_logistic_individual_IM(theta=theta_vec, Y=sum(Y), n=sum(n_i), M=50000)
lines(theta_vec, optimal_IM, col="black")

#############################################################################
## Compute Ryan's proposed naive IM
#############################################################################

# function for naive IM
pl.naive <- function(theta, theta_MLE, J_ni){
  s.J <- sum(J_ni)
  theta.n <- sum(J_ni*theta_MLE) / s.J
  1 - pchisq( ((theta.n -theta)^2)* s.J , 1)
}

# Plotting the naive IM
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni)
lines(theta_vec, naive_IM, col="red")

#############################################################################
## Compute the proposed middle ground IM with importance weights -- scalar theta, intercept only
#############################################################################

# work likelihood
work.likeli <- function(theta, theta_MLE, J_ni){
  if(any(J_ni == 0)) {
    if(sum(J_ni!=0) == 0) {
      ll <- NA
    } else {
      keep <- which(J_ni != 0)
      s.J <- sum(J_ni[keep])
      theta.n <- sum(theta_MLE[keep]*J_ni[keep]) / s.J
      ll <- exp(-(theta.n -theta)^2 * s.J / 2)
    }
    # ll <- NA
  } else {
    s.J <- sum(J_ni)
    theta.n <- sum(theta_MLE*J_ni) / s.J
    ll <- exp(-(theta.n -theta)^2 * s.J / 2)
  }
  return(ll)
}

# function for middle ground IM
pl.mid <- function(theta, theta.dag, theta_MLE, k, J_ni, n_i, M){
  Ysim_list <- lapply(1:length(theta.dag), function(t) sapply(1:k, function(i) {
    sapply(1:M, function(sim) rbinom(1, size=n_i[i], prob=expit(theta.dag[t])) )
    }) )
  
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    closest <- which.min(sum(abs(theta[l,] - theta.dag)))
    Ysim <- Ysim_list[[closest]]
    s_theta.l <- rowSums(sapply(1:k, function(i) {
      dbinom(Ysim[,i], size=n_i[i], prob=rep(expit(theta[l]), M), log=TRUE)
      }))
    s_theta.dag <- rowSums(sapply(1:k, function(i) {
      dbinom(Ysim[,i], size=n_i[i], prob=rep(expit(theta.dag[closest]), M), log=TRUE)
    }))
    weights <- exp(s_theta.l - s_theta.dag)
    res <- sapply(1:M, function(sim){
      theta_MLE_sim <- sapply(1:k, function(i) {
        if(Ysim[sim,i]==0) return(-Inf)
        else if (Ysim[sim,i]==n_i[i]) return(Inf)
        else return(log(Ysim[sim,i]/n_i[i]/(1-Ysim[sim,i]/n_i[i])))
      })
      J_ni_sim <- n_i * exp(theta_MLE_sim) / (1+exp(theta_MLE_sim))^2
      J_ni_sim[theta_MLE_sim %in% c(Inf, -Inf)] <- 0
      
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni)) * weights[sim]
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting middle ground IM
# theta.dag <- sum(theta_MLE*J_ni) / sum(J_ni)
# theta.dag <- theta_vec[seq(1,length(theta_vec),20)]
meta_estimate <- sum(theta_MLE*J_ni) / sum(J_ni)
theta.dag <- c(meta_estimate, meta_estimate+1)
middle_IM <- pl.mid(theta=theta_vec, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/logistic-contours-I.pdf", height=5, width=5)
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")
lines(theta_vec, optimal_IM, col="black")
lines(theta_vec, naive_IM, col="red")
lines(theta_vec, middle_IM, col="blue")
dev.off()

#############################################################################
### Comparison between middle ground IM and "exact" middle ground IM -- scalar theta, intercept only
#############################################################################

plot(theta_vec, middle_IM, col="blue", type="l",
     xlab=expression(theta), ylab=expression(pi(theta)))

# function for exact middle ground IM
pl.mid.exact <- function(theta, theta_MLE, k, J_ni, n_i, covariates, M){
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    Ysim <- sapply(1:k, function(i) {
      sapply(1:M, function(sim) rbinom(1, size=n_i[i], prob=expit(theta[l])) )
    })
    res <- sapply(1:M, function(sim) {
      theta_MLE_sim <- sapply(1:k, function(i) {
        if(Ysim[sim,i]==0) return(-Inf)
        else if (Ysim[sim,i]==n_i[i]) return(Inf)
        else return(log(Ysim[sim,i]/n_i[i]/(1-Ysim[sim,i]/n_i[i])))
      })
      J_ni_sim <- n_i * exp(theta_MLE_sim) / (1+exp(theta_MLE_sim))^2
      J_ni_sim[theta_MLE_sim %in% c(Inf, -Inf)] <- 0
      
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni))
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting exact middle ground IM
middle_IM_exact <- pl.mid.exact(theta=theta_vec, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=10000)
lines(theta_vec, middle_IM_exact, col="green")

#############################################################################
### Empirical CDF of the middle ground contour evaluated at theta_0
#############################################################################

nsim <- 10000
samples <- naive_samples <- rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  Y <- sapply(1:k, function(i) rbinom(n=1, size=n_i[i], prob=expit(theta)) )
  theta_MLE <- log(Y/n_i/(1-Y/n_i))
  if(any(Y==0) | any(Y==n_i) | any(is.na(theta_MLE))) {
    set.seed(sim+10000)
    Y <- sapply(1:k, function(i) rbinom(n=1, size=n_i[i], prob=expit(theta)) )
    theta_MLE <- log(Y/n_i/(1-Y/n_i))
    if(any(Y==0) | any(Y==n_i) | any(is.na(theta_MLE))) {
      set.seed(sim+20000)
      Y <- sapply(1:k, function(i) rbinom(n=1, size=n_i[i], prob=expit(theta)) )
      theta_MLE <- log(Y/n_i/(1-Y/n_i))
    }
  }
  
  J_ni <- n_i * exp(theta_MLE) / (1+exp(theta_MLE))^2
  meta_estimate <- sum(theta_MLE*J_ni) / sum(J_ni)
  theta.dag <- c(meta_estimate, meta_estimate+1)
  
  samples[sim] <- pl.mid(theta=theta, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)
  naive_samples[sim] <- pl.naive(theta=theta, theta_MLE=theta_MLE, J_ni=J_ni)
}

redo <- which(is.na(samples) | samples==0)
for(sim in redo){
  print(sim)
  set.seed(sim+30000)
  Y <- sapply(1:k, function(i) rbinom(n=1, size=n_i[i], prob=expit(theta)) )
  if(any(Y==0) | any(Y==n_i)) next

  theta_MLE <- log(Y/n_i/(1-Y/n_i))
  J_ni <- n_i * exp(theta_MLE) / (1+exp(theta_MLE))^2
  meta_estimate <- sum(theta_MLE*J_ni) / sum(J_ni)
  theta.dag <- c(meta_estimate, meta_estimate+1)

  samples[sim] <- pl.mid(theta=theta, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)
  naive_samples[sim] <- pl.naive(theta=theta, theta_MLE=theta_MLE, J_ni=J_ni)
}

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/ecdf_logistic-I.pdf", height=5, width=5)
plot(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(samples<=t)), 
     type="l", ylab="CDF", xlab=expression(alpha))
lines(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(naive_samples<=t)), col="red")
abline(0,1)
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
dev.off()

save.image("/Users/ehector/Dropbox/Projects/IMfusion/simulations/RData files/bernoulli-settingI-20240724.RData")
