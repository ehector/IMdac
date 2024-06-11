#############################################################################
## Individual and optimal plausibility contours
#############################################################################

individual_pl <- function(theta, Y, sigma) {
  2*(1-pnorm(q=abs(Y-theta)/sigma, mean=0, sd=1))
  # 1-pchisq((Y-theta)^2/sigma^2, df=1)
}
optimal_pl <- function(theta, Y, sigma, n) {
  Y_sigma <- sum(Y/sigma^2) / sum(1/sigma^2)
  1-abs(2*pnorm(q=Y_sigma-theta, mean=0, sd = sqrt(1/sum(1/sigma^2)) ) - 1)
}
recreate_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  1-abs(2*pnorm( (sum(1/sigma^2))^(-1/2) * sum(-sign(Y-theta)/sigma * qnorm(ipl/2)) ) - 1)
}

#############################################################################
## P-value combiners
#############################################################################

tippett_pl <- function(theta, Y, sigma, n) {
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  min(ipl)
}
fisher_pl <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  1-pchisq(-2*sum(log(ipl)), df=2*n)
}
pearson_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  pgamma(-sum(log(1-ipl)), shape=n, rate=1)
}
MG_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  x <- sum(log(ipl/(1-ipl)))
  1/(1+exp(-x))
}
edgington_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  x <- sum(ipl)
  sum(sapply(c(0:floor(x)), function(k) (-1)^k * choose(n,k) * (x-k)^n) )/ factorial(n)
}
stouffer_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  pnorm(sum(qnorm(ipl)))
}
prod_pl <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  prod(ipl)
}
prod_transf_pl <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  z <- prod(ipl)
  z * sum(sapply(1:(n-1), function(k) (-log(z))^k / factorial(k)))
}
kolmogorov_pl <- function(theta, Y, sigma, n, r){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  if(r == -Inf){ # precise
    a_r <- n
    return(min(1, a_r*min(ipl)))
  }
  if(r == -1){ # not precise
    a_r <- exp(1)*log(n)
    return(min(1, a_r * n/sum(1/ipl)))
  }
  if(r == 0){ # asymptotically precise
    a_r <- exp(1)
    return( min(1, a_r*prod(ipl)^(1/n) ))
  }
  if(r == Inf){ # precise
    a_r <- 1
    return(min(1, a_r*max(ipl)))
  }
  if(r < -1 & r != -Inf){ # asymptotically precise
    a_r <- r/(r+1)*n^(1+1/r)
    return(min(1, a_r * (sum(ipl^r)/n)^(1/r) ))
  }
  if(r > -1 & r <1/(n-1) & r != 0){ # asymptotically precise
    a_r <- (r+1)^(1/r)
    return(min(1, a_r * (sum(ipl^r)/n)^(1/r) ))
  }
  if(r >= 1/(n-1) & r <= n-1){ #precise
    a_r <- (r+1)^(1/r)
    return(min(1, a_r * (sum(ipl^r)/n)^(1/r) ))
  }
  if(r >= (n-1) & r != Inf){
    a_r <- n^(1/r)
    return(min(1, a_r * (sum(ipl^r)/n)^(1/r) ))
  }
}
kolmogorov_BA_pl <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  2*min( c(1/2, n*min(ipl), 2*mean(ipl) ))
}
kolmogorov_BG_pl <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  2*min( c(1/2, n*min(ipl), exp(1)*prod(ipl)^(1/n) ))
}

#############################################################################
## Possibility contour literature combiners
#############################################################################

hose_pl  <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  min(1-(1-ipl)^2)
}
HH_pl  <- function(theta, Y, sigma, n, weight){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  min(1, min(ipl/weight))
}
valid_HH_pl  <- function(theta, Y, sigma, n, weight){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  1 - prod(1- weight * min(1, min(ipl/weight)))
}

#############################################################################
## p-to-E-to-p calibrator combiners
#############################################################################

e_to_p_AR <- function(theta, Y, sigma, n){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  min(1, prod( sapply(ipl, function(p) {
    if(p == 1) return(0)
    else return(p * (-log(p))^2 / (1-p+p*log(p)) )
  }) ))
}
e_to_p_kappa <- function(theta, Y, sigma, n, weight){
  ipl <- sapply(1:n, function(i) individual_pl(theta, Y[i], sigma[i]) )
  min(1, prod( sapply(1:n, function(i) {
    if(ipl[i] == 1) return(0)
    else return( ipl[i]^(1-weight[i])/weight[i] )
  }) ))
}
