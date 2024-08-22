#############################################################################
## Generate data for Gaussian example, create individual and optimal contours
#############################################################################

set.seed(600)
theta <- 0
k <- 3
n_i <- c(5,10,15)
J_ni <- n_i/c(1,2,4)^2
theta_MLE <- sapply(1:k, function(i) rnorm(1, mean=theta, sd=1/sqrt(J_ni[i])))
#range of values for theta
theta_vec <- seq(-2,2,0.01)

# Function for the individual IM contours. 
gaussian_individual_IM <- function(theta, theta_MLE, J_ni){
  2*(1-pnorm(q=abs(theta_MLE-theta)*sqrt(J_ni)))
}

# Plotting the individual IM contours
individual_IMs <- lapply(1:k, function(i) gaussian_individual_IM(theta=theta_vec, theta_MLE=theta_MLE[i], J_ni=J_ni[i]) )

# Plotting the optimal IM contour
optimal_IM <- gaussian_individual_IM(theta=theta_vec, theta_MLE=sum(theta_MLE*J_ni)/sum(J_ni), J_ni=sum(J_ni))

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/simple-normal-opt-ind.pdf", height=5, width=5)
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")
lines(theta_vec, optimal_IM, col="black")
dev.off()

#############################################################################
## Compute Ryan's proposed naive IM
#############################################################################

pl.naive <- function(theta, theta_MLE, J_ni){
  s.J <- sum(J_ni)
  theta.n <- sum(J_ni*theta_MLE) / s.J
  1 - pchisq( ((theta.n -theta)^2)* s.J , 1)
}

# Plotting the naive IM
naive_IM <- pl.naive(theta_vec, theta_MLE=theta_MLE, J_ni=J_ni)
lines(theta_vec, naive_IM, col="red")

#############################################################################
## Compute the proposed middle ground IM with importance weights
#############################################################################

# work likelihood
work.likeli <- function(theta, theta_MLE, J_ni){
  s.J <- sum(J_ni)
  theta.n <- sum(theta_MLE*J_ni) / s.J
  return(exp(-(theta.n -theta)^2 * s.J / 2))
}

# function for middle ground IM
pl.mid <- function(theta, theta.dag, theta_MLE, k, J_ni, M){
  Ysim_list <- lapply(1:length(theta.dag), function(t) sapply(1:k, function(i) rnorm(M, mean=theta.dag[t], sd=1/sqrt(J_ni[i])) ))
  
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    closest <- which.min(abs(theta[l] - theta.dag))
    Ysim <- Ysim_list[[closest]]
    res <- sapply(1:M, function(sim){
      weights <-
        exp(sum(sapply(1:k, function(i) {
          dnorm(Ysim[sim,i], mean=theta[l], sd=1/sqrt(J_ni[i]), log=TRUE) - dnorm(Ysim[sim,i], mean=theta.dag[closest], sd=1/sqrt(J_ni[i]), log=TRUE)
        })))
      
      theta_MLE_sim <- Ysim[sim,]
      J_ni_sim <- J_ni
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni)) * weights
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting middle ground IM
# theta.dag <- sum(theta_MLE*J_ni) / sum(J_ni)
# theta.dag <- theta_vec[seq(1,length(theta_vec),4)]
# meta_estimate <- sum(theta_MLE*J_ni) / sum(J_ni)
# theta.dag <- c(meta_estimate-1, meta_estimate, meta_estimate+0.5, meta_estimate+1)
theta.dag <- seq(-2,2,0.05)
middle_IM <- pl.mid(theta=theta_vec, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, M=10000)

pdf("/Users/ehector/Dropbox/Projects/IMfusion/paper/figures_meta/gaussian-valid-meta-inference.pdf", height=5, width=5)
conf_level <- 0.8
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")
lines(theta_vec, middle_IM, col="blue")
abline(h=1-conf_level, lty=2, col="black")
abline(v=theta_vec[order(abs(middle_IM-1+conf_level), decreasing=F)[1]], lty=2, col="black")
abline(v=theta_vec[order(abs(middle_IM-1+conf_level), decreasing=F)[2]], lty=2, col="black")
dev.off()

#############################################################################
### Comparison between middle ground IM and "exact" middle ground IM
#############################################################################

plot(theta_vec, middle_IM, col="blue", type="l",
     xlab=expression(theta), ylab=expression(pi(theta)))

# function for exact middle ground IM
pl.mid.exact <- function(theta, theta.dag, theta_MLE, k, J_ni, M){
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    Ysim <- sapply(1:k, function(i) rnorm(M, mean=theta[l], sd=1/sqrt(J_ni[i])) )
    res <- sapply(1:M, function(sim) {
      as.numeric( work.likeli(theta=theta[l], theta_MLE=Ysim[sim,], J_ni=J_ni) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni))
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting exact middle ground IM
middle_IM_exact <- pl.mid.exact(theta=theta_vec, theta_MLE=theta_MLE, k=k, J_ni=J_ni, M=10000)
lines(theta_vec, middle_IM_exact, col="green")

