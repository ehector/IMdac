#############################################################################
## Generate data for Gamma example, create individual and optimal contours
#############################################################################

set.seed(600)
theta <- 0.5
k <- 3
n_i <- c(5,10,15)
Y <- sapply(1:k, function(i) rgamma(1, n_i[i], rate=theta) )
theta_MLE <- n_i/Y
J_ni <- n_i/theta_MLE^2

#range of values for theta
theta_vec <- seq(0,2,0.01)

# Function for the individual IM contours. 
gamma_individual_IM <- function(theta, Y, n){
  z <- (theta*Y)^n * exp(-theta*Y)
  expint::gammainc(n, -n*lamW::lambertWm1(-z^(1/n)/n))/factorial(n-1) -
    expint::gammainc(n, -n*lamW::lambertW0(-z^(1/n)/n))/factorial(n-1) + 1
}

# Plotting the individual IM contours
individual_IMs <- lapply(1:k, function(i) gamma_individual_IM(theta=theta_vec, Y=Y[i], n=n_i[i]) )
plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")

# Plotting the optimal IM contour
optimal_IM <- gamma_individual_IM(theta=theta_vec, Y=sum(Y), n=sum(n_i))
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
## Compute the proposed middle ground IM with importance weights
#############################################################################

# work likelihood
work.likeli <- function(theta, theta_MLE, J_ni){
  s.J <- sum(J_ni)
  theta.n <- sum(theta_MLE*J_ni) / s.J
  return(exp(-(theta.n -theta)^2 * s.J / 2))
}

# function for middle ground IM
pl.mid <- function(theta, theta.dag, theta_MLE, k, J_ni, n_i, M){
  Ysim_list <- lapply(1:length(theta.dag), function(t) sapply(1:k, function(i) rgamma(M, shape=n_i[i], rate=theta.dag[t])) )
  
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    closest <- which.min(abs(theta[l] - theta.dag))
    Ysim <- Ysim_list[[closest]]
    res <- sapply(1:M, function(sim){
      weights <-
        prod(sapply(1:k, function(i) {
          (theta[l]/theta.dag[closest])^n_i[i] * exp( (theta.dag[closest] - theta[l])*Ysim[sim,i])
          }))
      
      theta_MLE_sim <- n_i/Ysim[sim,]
      J_ni_sim <- n_i/theta_MLE_sim^2
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni)) * weights
    })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting middle ground IM
# theta.dag <- sum(theta_MLE*J_ni) / sum(J_ni)
theta.dag <- theta_vec[seq(1,length(theta_vec),20)]
middle_IM <- pl.mid(theta=theta_vec, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=10000)

plot(theta_vec, individual_IMs[[1]], type="l", col="grey", xlab=expression(theta), ylab=expression(pi(theta)))
lines(theta_vec, individual_IMs[[2]],col="grey")
lines(theta_vec, individual_IMs[[3]],col="grey")
lines(theta_vec, optimal_IM, col="black")
lines(theta_vec, middle_IM, col="blue")

#############################################################################
### Comparison between middle ground IM and "exact" middle ground IM
#############################################################################

plot(theta_vec, middle_IM, col="blue", type="l",
     xlab=expression(theta), ylab=expression(pi(theta)))

# function for exact middle ground IM
pl.mid.exact <- function(theta, theta.dag, theta_MLE, k, J_ni, n_i, M){
  pl <- vector(length=length(theta))
  for(l in 1:length(theta)){
    Ysim <- sapply(1:k, function(i) rgamma(M, shape=n_i[i], rate=theta[l]))
    res <- sapply(1:M, function(sim) {
      theta_MLE_sim <- n_i/Ysim[sim,]
      J_ni_sim <- n_i/theta_MLE_sim^2
      as.numeric( work.likeli(theta=theta[l], theta_MLE=theta_MLE_sim, J_ni=J_ni_sim) <= 
                    work.likeli(theta=theta[l], theta_MLE=theta_MLE, J_ni=J_ni))
      })
    pl[l] <- mean(res)
  }
  return(pl)
}

# Plotting exact middle ground IM
middle_IM_exact <- pl.mid.exact(theta=theta_vec, theta.dag=theta.dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=10000)
lines(theta_vec, middle_IM_exact, col="green")

