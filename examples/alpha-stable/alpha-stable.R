## Conda set up to run Bayesflow (some of this from https://bayesflow.org/installation.html)
## build conda environment using
# ~$ conda create --prefix /usr/local/usrapps/$GROUP/$USER/env_bf -y python=3.10
## Then activate the environment and install bayesflow:
# ~$ conda activate /usr/local/usrapps/$GROUP/$USER/env_bf
# ~$ pip install bayesflow

## Load necessary packages, including IMdac
library(IMdac)
library(numDeriv)
library(reticulate)
library(stabledist)

## Set-up python directory
python_dir <- "/usr/local/usrapps/$GROUP/$USER/env_bf" # for Hazel
use_condaenv(python_dir)

## The path where the amortizer lives on your computer
amortizer_path <- "IMdac/examples/alpha-stable/trained_alpha_stable_amortizer"

## Import numpy and bayesflow
np <- import("numpy")
bf <- import("bayesflow")

## The alpha-stable negative log-likelihood. Takes as input a vector "par"
## consisting of beta, c, and mu
alpha_stable_nll <- function(par, alpha, y){
  beta <- par[1]
  c <- par[2]
  mu <- par[3]
  if(beta < -1 | beta > 1) return(-Inf)
  else {
    return( -sum(sapply(1:length(y), function(i) dstable(y[i], alpha=alpha, beta=beta, gamma=c, delta=mu, pm=1, log=T) )) ) 
  }
}

## Generate the data
sim <- 600
set.seed(sim)
n_i <- rep(500,4)
k <- length(n_i)
alpha <- 0.5
beta <- 0
c <- 0.5
mu <- 0
theta <- c(beta, c, mu)
p <- length(theta)
Y <- lapply(1:k, function(i) {
  U <- runif(n_i[i], -pi/2, pi/2)
  W <- rexp(n_i[i], 1)
  zeta <- -beta * tan(pi*alpha/2)
  if(alpha==1){
    xi <- pi/2
    X <- 1/xi*( ( pi/2 + beta*U)*tan(U) - beta*log( (pi/2*W*cos(U))/(pi/2+beta*U) ))
    Y <- c*X + 2*beta*c*log(c)/pi + mu
  } else{
    xi <- atan(-zeta)/alpha
    X <- (1+zeta^2)^(1/(2*alpha)) * 
      (sin(alpha*(U+zeta)))/(cos(U))^(1/alpha) * 
      ( (cos(U - alpha*(U+xi)))/W)^((1-alpha)/alpha)
    Y <- c*X + mu
  }
  return(Y)
  } )
rm(beta, c, mu)

## Dimension of summary statistic in amortizer
num_summary_dim <- 10L

## Train the amortizer. Skip this step if using pre-trained amortizer
amortizer <- train_alpha_stable_emulator(n=n_i[1], num_summary_dim=num_summary_dim, alpha=alpha)
amortizer$save_weights(amortizer_path)

## If the amortizer is trained, load the trained amortizer
amortizer <- load_trained_alpha_stable_emulator(n=n_i[1], num_summary_dim=num_summary_dim, alpha=alpha, path=amortizer_path)

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
        large_n_estimate[3]-3*sqrt(diag(solve(s.J)))[3])

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
  alpha_stable_pl_mid_all_profile_NN(
    theta=theta_q_vec, alpha=alpha, theta_dag=large_n_estimate,
    large_n_estimate=large_n_estimate, k, s_J=s.J, n_i=n_i, p=p,
    M=3000, n_samples=n_samples, amortizer=amortizer)
time <- proc.time() - time

## the full data MLE and observed Fisher information
time_full <- proc.time()
MLE_full_optim <- optim(par=rep(0.5,p), fn=alpha_stable_nll, alpha=alpha, y=unlist(Y), method="Nelder-Mead", control=list(maxit=5000))
MLE_full <- MLE_full_optim$par
J_ni_full <- numDeriv::hessian(function(x) alpha_stable_nll(x, unlist(Y), alpha=alpha), 
                               x=MLE_full)
time_full <- proc.time() - time_full
