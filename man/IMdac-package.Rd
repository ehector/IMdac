\name{IMdac-package}
\alias{IMdac-package}
\alias{IMdac}
\docType{package}
\title{\packageTitle{IMdac}}
\description{\packageDescription{IMdac}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{IMdac}
  \packageIndices{IMdac}

  This is a repository for the R package to perform possibilistic divide-and-conquer analysis in the inferential models framework. The R package's main file is the R/dac_funcs-20250221.R file. It defines the function for the divide-and-conquer possibility contour and helper functions. The examples folder contains the code to reproduce the alpha-stable and g-and-k examples. Due to the use of Bayesflow, some set-up of Python is needed that is system specific. Instructions are available in the comments at the top of the two scripts: "alpha-stable.R" and "g-and-k.R".
}
\author{
  \packageAuthor{IMdac}
  
  Maintainer: \packageMaintainer{IMdac}
}
\references{
  # Installation

The IMdac R package can be installed in one of two ways:

- from the downloaded gzipped tarball as R CMD INSTALL IMdac_1.0-1.tar.gz

- from the downloaded and renamed IMdac folder as R CMD build IMdac and R CMD INSTALL IMdac_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the IMdac R package, please consider citing the relevant manuscript: E.C. Hector, L. Cella and R. Martin (2025+). Divide-and-conquer with finite sample sizes: valid and efficient possibilistic inference.

# References

Hector, E. C., Tang, L., Zhou, L., and Song, P. X.-K. (2024). Handbook on Bayesian, Fiducial and Frequentist Inference, chapter Data integration and fusion in the Bayesian and Frequentist frameworks. Chapman and Hall/CRC Press.

R. Martin and C. Liu. (2015). Inferential models: reasoning with uncertainty. CRC Press.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
}
\examples{
###############################################################################
## Exponential example
###############################################################################

set.seed(600)
## the true value of the uncertain Theta
theta <- 0.5

## the number of blocks
k <- 3

## the size of each block
n_i <- c(5,10,15)

## the data in each block
Y <- sapply(1:k, function(i) rgamma(1, n_i[i], rate=theta) )

## the MLE in each block
theta_MLE <- n_i/Y

## the observed Fisher information in each block
J_ni <- n_i/theta_MLE^2

## range of values for theta
theta_vec <- seq(0,2,0.01)

## compute the large-n estimator
s.J <- sum(J_ni)
large_n_estimate <- sum(J_ni*theta_MLE) / s.J

## Function for the exact IM contours
gamma_individual_IM <- function(theta, Y, n){
  z <- (theta*Y)^n * exp(-theta*Y)
  expint::gammainc(n, -n*lamW::lambertWm1(-z^(1/n)/n))/factorial(n-1) -
    expint::gammainc(n, -n*lamW::lambertW0(-z^(1/n)/n))/factorial(n-1) + 1
}

## the individual IM contours
individual_IMs <- lapply(1:k, function(i) gamma_individual_IM(theta=theta_vec, Y=Y[i], n=n_i[i]) )

## the optimal IM contour
optimal_IM <- gamma_individual_IM(theta=theta_vec, Y=sum(Y), n=sum(n_i))

## the large-n contour
naive_IM <- sapply(theta_vec, function(theta) {
  1 - pchisq( ((large_n_estimate -theta)^2)* s.J , 1)
  })

## the valid divide-and-conquer contour
theta_dag <- matrix(meta_estimate, nrow=1)
middle_IM <- exponential_pl_mid_P(theta=theta_vec, theta_dag=theta_dag, theta_MLE=theta_MLE, k=k, J_ni=J_ni, n_i, M=50000)

###############################################################################
## Gaussian example
###############################################################################

set.seed(600)

## the true value of the uncertain Theta
theta <- 0

## the number of blocks
k <- 3

## the size of each block
n_i <- c(5,10,15)

## the MLE and observed Fisher information in each block
J_ni <- n_i/c(2,2,2)^2
theta_MLE <- sapply(1:k, function(i) rnorm(1, mean=theta, sd=1/sqrt(J_ni[i])))

## range of values for theta
theta_vec <- seq(-2,2,0.01)

## Function for the exact IM contours
gaussian_individual_IM <- function(theta, theta_MLE, J_ni){
  2*(1-pnorm(q=abs(theta_MLE-theta)*sqrt(J_ni)))
}

## the individual IM contours
individual_IMs <- lapply(1:k, function(i) gaussian_individual_IM(theta=theta_vec, theta_MLE=theta_MLE[i], J_ni=J_ni[i]) )

## the optimal IM contour
optimal_IM <- gaussian_individual_IM(theta=theta_vec, theta_MLE=sum(theta_MLE*J_ni)/sum(J_ni), J_ni=sum(J_ni))

## the large-n contour
naive_IM <- sapply(theta_vec, function(theta) {
  1 - pchisq( ((large_n_estimate -theta)^2)* s.J , 1)
  })

## the valid divide-and-conquer contour
middle_IM <- naive_IM
}
