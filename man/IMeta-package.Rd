\name{IMeta-package}
\alias{IMeta-package}
\alias{IMeta}
\docType{package}
\title{\packageTitle{IMeta}}
\description{\packageDescription{IMeta}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{IMeta}
  \packageIndices{IMeta}

  This is a repository for the R package to perform possibilistic meta-analysis in the inferential models framework. The R package's main files are:

- src/compute_funcs-20240828.cpp: this file defines the Rcpp functions that compute the possibility contours for the individual, optimal, meta-analysis and equilibrated possibility contours.

- R/meta_analysis_funcs-20240902.R: this file defines the function for the meta-analysis possibility contour and helper functions.
}
\author{
  \packageAuthor{IMeta}
  
  Maintainer: \packageMaintainer{IMeta}
}
\references{
  # Installation

The IMeta R package can be installed in one of two ways:

- from the downloaded gzipped tarball as R CMD INSTALL IMeta_1.0-1.tar.gz

- from the downloaded and renamed IMeta folder as R CMD build IMeta and R CMD INSTALL IMeta_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the IMeta R package, please consider citing the relevant manuscript: E.C. Hector, L. Cella and R. Martin (2024+). Efficient possibilistic meta-analysis with finite-sample validity guarantees.

# References

G.V. Glass. (1976). Primary, secondary, and meta-analysis of research. Educational Researcher, 5(10):3–8.

R. Martin and C. Liu. (2015). Inferential models: reasoning with uncertainty. CRC Press.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
}
\examples{
########### Lognormal example ########### 

## Set parameters
n_i <- c(12,14,10,11,13)
k <- length(n_i)
mu <- 1
gamma <- 1.5
theta <- c(mu, gamma)
p <- length(theta)

## Generate the data
set.seed(600)
Y <- lapply(1:k, function(i) rlnorm(n=n_i[i], meanlog=mu, sdlog=sqrt(gamma)) )

## Compute the MLE's and the observed Fisher information
theta_MLE <- t(sapply(1:k, function(i) lognormal_MLE(Y[[i]]) ))
J_ni <- lapply(1:k, function(i) lognormal_FI(Y[[i]], n_i[i], theta_MLE[i,]))

## Set the range of values for theta
s.J <- Reduce("+",J_ni)
meta_estimate <- drop(solve(s.J) \%*\% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]\%*\%theta_MLE[i,]) ))
theta_vec <- as.matrix(expand.grid(
  seq(meta_estimate[1]-5*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+5*sqrt(diag(solve(s.J)))[1],length=400), 
  seq(meta_estimate[2]-5*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+5*sqrt(diag(solve(s.J)))[2],length=400)) )
  
## Compute the meta-analysis possibility contour
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni, p=p)

## Compute the equilibrated possibility contour
theta_dag <- as.matrix(expand.grid(
  c(meta_estimate[1] - 1*sqrt(diag(solve(s.J)))[1], meta_estimate[1], meta_estimate[1] + 1*sqrt(diag(solve(s.J)))[1]),
  c(meta_estimate[2] - 1*sqrt(diag(solve(s.J)))[2], meta_estimate[2], meta_estimate[2] + 1*sqrt(diag(solve(s.J)))[2])
))
middle_IM <- lognormal_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i=n_i, M=50000)

########### Logistic regression example ########### 

## Set parameters
theta <- 0.5
k <- 3
n_i <- c(50,75,100)
theta <- c(-0.5,0.3)
p <- length(theta)

## Generate the data
set.seed(600)
covariates <- lapply(1:k, function(i) cbind(1, rnorm(n_i[i],0,1)) )
Y <- sapply(1:k, function(i) rbinom(n=n_i[i], size=1, prob=expit(covariates[[i]]\%*\%theta)) )

## Compute the MLE's and the observed Fisher information
theta_MLE <- matrix(0,k,p)
J_ni <- list()
for(i in 1:k){
  glm_i <- glm(Y[[i]] ~ 0 + covariates[[i]], family="binomial")
  theta_MLE[i,] <- as.vector(coef(glm_i))
  mu_i <- expit(covariates[[i]]%*%theta_MLE[i,])
  J_ni[[i]] <- Reduce("+", lapply(1:n_i[i], function(j) covariates[[i]][j,] \%o\% covariates[[i]][j,] * mu_i[j]*(1-mu_i[j])))
}

## Set the range of values for theta
s.J <- Reduce("+",J_ni)
meta_estimate <- drop(solve(s.J) \%*\% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]\%*\%theta_MLE[i,]) ))
theta_vec <- as.matrix(expand.grid(
  seq(meta_estimate[1]-2*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+2*sqrt(diag(solve(s.J)))[1],length=150), 
  seq(meta_estimate[2]-2*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+2*sqrt(diag(solve(s.J)))[2],length=150)) )

## Compute the meta-analysis possibility contour
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni, p=p)

## Compute the equilibrated possibility contour
theta_dag <- matrix(meta_estimate, nrow=1)
middle_IM <- logistic_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i,
                             covariates=covariates, M=50000)

########### Count regression example ########### 

## Set parameters
k <- 3
n_i <- c(15,20,25)
theta <- c(-0.2,0.5,0.8)
p <- length(theta)

## Generate the data
set.seed(600)
covariates <- lapply(1:k, function(i) cbind(1, rnorm(n_i[i],0,1), rnorm(n_i[i],0,1)) )
Y <- lapply(1:k, function(i) rpois(n=n_i[i], lambda=exp(covariates[[i]]\%*\%theta)) )

## Compute the MLE's and the observed Fisher information
theta_MLE <- matrix(0,k,p)
J_ni <- list()
for(i in 1:k){
  glm_i <- glm(Y[[i]] ~ 0 + covariates[[i]], family="poisson")
  theta_MLE[i,] <- as.vector(coef(glm_i))
  mu_i <- exp(covariates[[i]]\%*\%theta_MLE[i,])
  J_ni[[i]] <- Reduce("+", lapply(1:n_i[i], function(j) covariates[[i]][j,] \%o\% covariates[[i]][j,] * mu_i[j]))
}

## Set the range of values for theta
s.J <- Reduce("+",J_ni)
meta_estimate <- drop(solve(s.J) \%*\% Reduce("+", lapply(1:length(J_ni), function(i) J_ni[[i]]\%*\%theta_MLE[i,]) ))
theta_vec <- as.matrix(expand.grid(
  seq(meta_estimate[1]-4*sqrt(diag(solve(s.J)))[1], meta_estimate[1]+4*sqrt(diag(solve(s.J)))[1],length=100), 
  seq(meta_estimate[2]-4*sqrt(diag(solve(s.J)))[2], meta_estimate[2]+4*sqrt(diag(solve(s.J)))[2],length=100),
  seq(meta_estimate[3]-4*sqrt(diag(solve(s.J)))[3], meta_estimate[3]+4*sqrt(diag(solve(s.J)))[3],length=100)) )

## Compute the meta-analysis possibility contour
naive_IM <- pl.naive(theta=theta_vec, theta_MLE=theta_MLE, J_ni=J_ni, p=p)

## Compute the equilibrated possibility contour
theta_dag <- matrix(meta_estimate, nrow=1)

middle_IM <- poisson_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, covariates=covariates, M=50000)


}
