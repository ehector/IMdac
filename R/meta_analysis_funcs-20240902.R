expit <- function(z) exp(z)/(1+exp(z))

pl.naive <- function(theta, theta_MLE, J_ni, p){
  s.J <- Reduce("+",J_ni)
  theta.n <- solve(s.J) %*% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]%*%theta_MLE[i,]) )
  apply(theta, 1, function(t) drop(1 - pchisq( t(theta.n -t) %*% s.J %*% (theta.n -t), p)))
}