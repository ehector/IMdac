## Conda set up to run Bayesflow (some of this from https://bayesflow.org/installation.html)
## build conda environment using
# ~$ conda create --prefix /usr/local/usrapps/$GROUP/$USER/env_bf -y python=3.10
## Then activate the environment and install bayesflow:
# ~$ conda activate /usr/local/usrapps/$GROUP/$USER/env_bf
# ~$ pip install bayesflow

## Load necessary packages, including IMdac
library(gk)
library(IMdac)
library(numDeriv)
library(reticulate)
Rcpp::sourceCpp("IMdac/examples/g-and-k/g-and-k-functions.cpp")

## Set-up python directory
python_dir <- "/usr/local/usrapps/$GROUP/$USER/env_bf" # for Hazel
use_condaenv(python_dir)

## The path where the amortizer lives on your computer
amortizer_path <- "IMdac/examples/g-and-k/trained_g_and_k_amortizer"

## Import numpy and bayesflow
np <- import("numpy")
bf <- import("bayesflow")

## The g-and-k negative log-likelihood. 
## Some of this is borrowed from the "dgk" function in the R package "gk".
dgk <- function(x, A, B, g, k, c) {
  z = pgk(x, A, B, g, k, c, zscale = TRUE)
  
  return(stats::dnorm(z, log = TRUE) - Cpp_Qgk_log_deriv(z, A, B, g, k, c))
}

## Some of this is borrowed from the "z2gk" function in the R package "gk".
z2gk <- function (z, A, B, g, k, c) {
  n = max(length(z), length(A), length(B), length(g), length(k), 
          length(c))
  zeros = rep(0, n)
  z = z + zeros
  A = A + zeros
  B = B + zeros
  g = g + zeros
  k = k + zeros
  c = c + zeros
  z_squared = z^2
  term1 = (1 + c * tanh(g * z/2))
  term2 = z * (1 + z_squared)^k
  term1[g == 0] = 1
  zbig = which(is.infinite(z_squared))
  term2[zbig] = sign(z[zbig]) * abs(z[zbig])^(1 + 2 * k[zbig])
  return(A + B * term1 * term2)
}

## Some of this is borrowed from the "rgk" function in the R package "gk".
rgk <- function (n, A, B, g, k, c) {
  z2gk(stats::rnorm(n), A, B, g, k, c)
}

## Some of this is borrowed from the "pgk_scalar" function in the R package "gk".
pgk_scalar <- function (q, A, B, g, k, c, zscale) {
  toroot = function(p) {
    z2gk(p, A, B, g, k, c) - q
  }
  z = tryCatch(stats::uniroot(toroot, interval = c(-5, 5), 
                              extendInt = "upX", check.conv = TRUE)$root, error = function(cond) {
                                Inf * (q - A)
                              })
  if (zscale) {
    return(z)
  }
  return(stats::pnorm(z))
}

## Some of this is borrowed from the "pgk" function in the R package "gk".
pgk <- function (q, A, B, g, k, c, zscale) {
  sapply(q, function(x) pgk_scalar(x, A, B, g, k, c, zscale))
}

## Generate the data
sim <- 600
set.seed(sim)
n_i <- rep(50,4)
k <- length(n_i)
A <- 3
B <- 1
g <- 2
kk <- 0.5
theta <- c(A, B, g, kk)
rm(A, B, g, kk)
p <- length(theta)
Y <- lapply(1:k, function(i) rgk(n=n_i[i], A=theta[1], B=theta[2], g=theta[3], k=theta[4], c=0.8) )

## Dimension of summary statistic in amortizer
num_summary_dim <- 10L

## Train the amortizer. Skip this step if using pre-trained amortizer
amortizer <- train_g_and_k_emulator(n=n_i[1], num_summary_dim=num_summary_dim)
amortizer$save_weights(amortizer_path)

## If the amortizer is trained, load the trained amortizer
amortizer <- load_trained_g_and_k_emulator(n=n_i[1], num_summary_dim=num_summary_dim, path=amortizer_path)

## the MLE and observed Fisher information in each block
n_samples <- 10000L
emulator_samples <- lapply(1:k, function(i) {
  ## do some set up to pass the data into the amortizer
  conf_data <- list()
  conf_data$summary_conditions <- array(Y[[i]], dim=c(1,n_i[i],1))
  conf_data$direct_conditions <- NULL
  emulator_samples <- amortizer$sample(conf_data, n_samples=n_samples)
  
  ## compute the MLE and variance as the mean and variance of draws from the emulator
  MLE <- colMeans(emulator_samples)
  var <- var(emulator_samples)
  return(list(MLE=MLE, var=var))
})
theta_MLE <- t(sapply(1:k, function(i) {
  emulator_samples[[i]]$MLE
  }))
J_ni <- lapply(1:k, function(i) {
  solve(emulator_samples[[i]]$var)
})

## compute the large-n estimator
s.J <- Reduce("+",J_ni)
large_n_estimate <- drop(solve(s.J) %*% Reduce("+", lapply(1:k, function(i) J_ni[[i]]%*%theta_MLE[i,]) ))

## lower bound for range of values for theta
lb <- c(large_n_estimate[1]-3*sqrt(diag(solve(s.J)))[1], max(1e-5, large_n_estimate[2]-3*sqrt(diag(solve(s.J)))[2]), 
        large_n_estimate[3]-3*sqrt(diag(solve(s.J)))[3], max(1e-5, large_n_estimate[4]-3*sqrt(diag(solve(s.J)))[4]))

## number of theta values in the sequence
theta_length <- 10000

naive_IM <- theta_q_vec <- matrix(0, theta_length, p)
time_naive <- list()
for(q in 1:p){
  ## for each parameter, set up the theta sequence
  print(paste0("q=",q))
  theta_q_vec[,q] <- seq(lb[q],large_n_estimate[q]+3*sqrt(diag(solve(s.J)))[q], length=theta_length)
  
  ## for each parameter, compute the large-n contour
  time_naive[[q]] <- proc.time()
  naive_IM[,q] <- sapply(theta_q_vec[,q], function(t) drop(1 - pchisq( s.J[q,q]*(large_n_estimate[q]-t)^2, 1)) )
  time_naive[[q]] <- proc.time() - time_naive[[q]]
}

## the valid divide-and-conquer contour
time <- proc.time()
middle_IM <-
  g_and_k_pl_mid_all_profile_NN(theta=theta_q_vec, theta_dag=large_n_estimate,
                                large_n_estimate=large_n_estimate, k, s_J=s.J, n_i=n_i, p=p,
                                M=3000, n_samples=n_samples, amortizer=amortizer)
time <- proc.time() - time

## set-up to compute the full data MLE. 
## Some of this is borrowed from the "fdsa" function in the R package "gk".
n_it <- 1000
theta_min <- c(-5,1e-5,-5,1e-5)
theta_max <- c(5,5,5,5)
theta0 <- c(0,1,0,1e-5)
density_sample <- dgk(unlist(Y), theta0[1], theta0[2], theta0[3], theta0[4], 0.8)
c0 <- sd(density_sample)/sqrt(sum(n_i))
c0 <- pmin(c0, (theta_max - theta_min)/2)

## compute the full data MLE and observed Fisher information
time_full <- proc.time()
MLE_full <- g_and_k_MLE(unlist(Y), 1, 100, 1.0, c0, 0.49, p, theta0=theta0, 
                        theta_min=theta_min, theta_max=theta_max, maxit=n_it, 
                        tol=1e-6)
J_ni_full <- numDeriv::hessian(function(x) {
  sum(-dgk(unlist(Y), A=x[1], B=x[2], g=x[3], k=x[4], c=0.8))}, 
  x=MLE_full)
time_full <- proc.time() - time_full

