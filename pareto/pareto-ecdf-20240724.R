#############################################################################
## Generate data for Pareto example, create individual and optimal contours -- scalar theta, intercept only
#############################################################################

set.seed(600)
k <- 3
n_i <- c(10,15,20)
theta <- 10
pareto_QF <- function(u, theta) (1-u)^(-1/theta)
Y <- sapply(1:k, function(i) pareto_QF(runif(n_i[i], 0, 1), theta=theta) )
theta_MLE <- sapply(1:k, function(i) n_i[i]/sum(log(Y[[i]])))
J_ni <- n_i / theta_MLE^2

#range of values for theta
meta_estimate <- sum(theta_MLE*J_ni) / sum(J_ni)
theta_vec <- seq(meta_estimate-2/sum(J_ni), meta_estimate+2/sum(J_ni), length=5000)

pareto_individual_IM <- function(theta, Y, n){
  exp(n*log(theta)+((n/sum(log(Y))+1)-(theta+1))*sum(log(Y)) - n*log(n/sum(log(Y))) )
}
valid_pareto_individual_IM <- function(theta, Y, n, M){
  observed_contour <- sapply(theta, function(t) pareto_individual_IM(theta=t, Y=Y, n=n))
  pl <- rep(0, length(theta))
  for(t in 1:length(theta)){
    Ysim <- t(sapply(1:M, function(sim) pareto_QF(runif(n, 0, 1), theta=theta[t]) ))
    res <- sapply(1:M, function(sim) {
      as.numeric( pareto_individual_IM(theta=theta[t], Y=Ysim[sim,], n=n) <= observed_contour[t])
    })
    pl[t] <-   mean(res)
  }
  return(pl)
}

# Individual IM contours
individual_IMs <- lapply(1:k, function(i) valid_pareto_individual_IM(theta=theta_vec, Y=Y[[i]], n=n_i[i], M=50000) )

# Optimal IM contour
optimal_IM <- valid_pareto_individual_IM(theta=theta_vec, Y=unlist(Y), n=sum(n_i), M=50000)

#############################################################################
## Ryan's proposed naive IM
#############################################################################

# function for naive IM
pl.naive <- function(theta, theta_MLE, J_ni){
  s.J <- sum(J_ni)
  theta.n <- sum(J_ni*theta_MLE) / s.J
  1 - pchisq( ((theta.n -theta)^2)* s.J , 1)
}

# Naive IM
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni)

#############################################################################
## Proposed middle ground IM with importance weights -- scalar theta, intercept only
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
  Ysim_list <- lapply(1:length(theta.dag), function(t) lapply(1:k, function(i) {
    t(sapply(1:M, function(sim) pareto_QF(runif(n_i[i], 0, 1), theta=theta.dag[t]) ))
    }) )
  
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    closest <- which.min(sum(abs(theta[l] - theta.dag)))
    Ysim <- Ysim_list[[closest]]
    s_theta.l <- rowSums(sapply(1:k, function(i) {
      n_i[i]*log(theta[l]) - (theta[l]+1)*rowSums(log(Ysim[[i]]))
      }))
    s_theta.dag <- rowSums(sapply(1:k, function(i) {
      n_i[i]*log(theta.dag[closest]) - (theta.dag[closest]+1)*rowSums(log(Ysim[[i]]))
    }))
    weights <- exp(s_theta.l - s_theta.dag)
    res <- sapply(1:M, function(sim){
      theta_MLE_sim <- sapply(1:k, function(i) {
        if(any(Ysim[[i]][sim,]==0)) return(-Inf)
        else return( n_i[i]/sum(log(Ysim[[i]][sim,])) )
      })
      J_ni_sim <-  n_i / theta_MLE_sim^2
      J_ni_sim[theta_MLE_sim == -Inf] <- 0
      
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni)) * weights[sim]
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Middle ground IM
# theta.dag <- c(meta_estimate-2, meta_estimate-1, meta_estimate, meta_estimate+1, meta_estimate+2)
# middle_IM <- pl.mid(theta=theta_vec, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)

# function for middle ground IM
pl.mid.exact <- function(theta, theta_MLE, k, J_ni, n_i, M){
  
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    
    Ysim <- lapply(1:k, function(i) {
      t(sapply(1:M, function(sim) pareto_QF(runif(n_i[i], 0, 1), theta=theta[l]) ))
    })
    
    res <- sapply(1:M, function(sim){
      theta_MLE_sim <- sapply(1:k, function(i) {
        if(any(Ysim[[i]][sim,]==0)) return(-Inf)
        else return( n_i[i]/sum(log(Ysim[[i]][sim,])) )
      })
      J_ni_sim <-  n_i / theta_MLE_sim^2
      J_ni_sim[theta_MLE_sim == -Inf] <- 0
      
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <=
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni))
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

middle_IM_exact <- pl.mid.exact(theta=theta_vec, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/pareto-contours.pdf", height=5, width=5)
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)), ylim=c(0,1))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")
lines(theta_vec, optimal_IM, col="black")
lines(theta_vec, naive_IM, col="red")
lines(theta_vec, middle_IM_exact, col="blue")
dev.off()

#############################################################################
### Empirical CDF of the middle ground contour evaluated at theta_0
#############################################################################

nsim <- 10000
samples <- naive_samples <- rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  Y.sim <- sapply(1:k, function(i) pareto_QF(runif(n_i[i], 0, 1), theta=theta) )
  theta_MLE.sim <- sapply(1:k, function(i) n_i[i]/sum(log(Y.sim[[i]])))
  J_ni.sim <- n_i / theta_MLE.sim^2

  samples[sim] <- pl.mid.exact(theta=theta, theta_MLE=theta_MLE.sim, k=k, J_ni=J_ni.sim, n_i, M=50000)
  naive_samples[sim] <- pl.naive(theta=theta, theta_MLE=theta_MLE.sim, J_ni=J_ni.sim)
}

save.image("pareto-ecdf.RData")

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/ecdf_pareto.pdf", height=5, width=5)
plot(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(samples<=t)),
     type="l", ylab="CDF", xlab=expression(alpha))
lines(seq(0,1,length=100), sapply(seq(0,1,length=100), function(t) mean(naive_samples<=t)), col="red")
abline(0,1)
abline(h=0, lty=2, col="gray")
abline(h=1, lty=2, col="gray")
dev.off()
