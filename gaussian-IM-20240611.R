#############################################################################
## Call R packages, source functions
#############################################################################

library(ggplot2)
source("gaussian-combiners-funcs-20240611.R")

#############################################################################
## Setup
#############################################################################

sigma <- c(1,2,4)
n <- length(sigma)
theta <- 0

set.seed(20240401)
Y <- sapply(1:n, function(i) rnorm(1, mean=theta, sd=sigma[i]))
theta_vec <- seq(floor(min(Y)-3), ceiling(max(Y)+3), length=10000)

#############################################################################
## Individual and optimal plausibility contours
#############################################################################

eval_individual1 <- sapply(theta_vec, function(t) individual_pl(t, Y[1], sigma[1]))
eval_individual2 <- sapply(theta_vec, function(t) individual_pl(t, Y[2], sigma[2]))
eval_individual3 <- sapply(theta_vec, function(t) individual_pl(t, Y[3], sigma[3]))
eval_optimal <- sapply(theta_vec, function(t) optimal_pl(t, Y, sigma, n))

#############################################################################
## P-value combiners
#############################################################################

eval_tippett <- sapply(theta_vec, function(t) tippett_pl(t, Y, sigma, n)) / 
  max(sapply(theta_vec, function(t) tippett_pl(t, Y, sigma, n)))
eval_fisher <- sapply(theta_vec, function(t) fisher_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) fisher_pl(t, Y, sigma, n)))
eval_pearson <- sapply(theta_vec, function(t) pearson_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) pearson_pl(t, Y, sigma, n)))
eval_MG <- sapply(theta_vec, function(t) MG_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) MG_pl(t, Y, sigma, n)))
eval_edgington <- sapply(theta_vec, function(t) edgington_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) edgington_pl(t, Y, sigma, n)))
eval_stouffer <- sapply(theta_vec, function(t) stouffer_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) stouffer_pl(t, Y, sigma, n)))
eval_prod <- sapply(theta_vec, function(t) prod_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) prod_pl(t, Y, sigma, n)))
eval_prod_transf <- sapply(theta_vec, function(t) prod_transf_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) prod_transf_pl(t, Y, sigma, n)))

eval_kolmogorov_min <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-Inf))/
  max(sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-Inf)))
eval_kolmogorov_max <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=Inf))
eval_kolmogorov_harmonic <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-1))/
  max(sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-1)))
eval_kolmogorov_geometric <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=0))/
  max(sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=0)))
eval_kolmogorov_arithmetic <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=1))
eval_kolmogorov_neg <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-1.5))
# eval_kolmogorov_n <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=-0.5))
# eval_kolmogorov_p <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=1.5))
eval_kolmogorov_pos <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=3))

eval_kolmogorov_BA <- sapply(theta_vec, function(t) kolmogorov_BA_pl(t, Y, sigma, n))
eval_kolmogorov_BG <- sapply(theta_vec, function(t) kolmogorov_BG_pl(t, Y, sigma, n))

#############################################################################
## Optimization combiners
#############################################################################

eval_median <- sapply(theta_vec, function(t) median_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) median_pl(t, Y, sigma, n)))
eval_p25 <- sapply(theta_vec, function(t) p25_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) p25_pl(t, Y, sigma, n)))

int <- sum(2*sigma/sqrt(pi))
opt <- optim(par=0, fn=kld_opt, Y=Y, sigma=sigma, n=n, int=int, ll=-10, ul=10, method="Brent", lower=-500, upper=500)
eval_kolmogorov_opt <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=opt$par))

opt <- optim(par=0, fn=tv_opt, Y=Y, sigma=sigma, n=n, theta=theta_vec, method="Brent", lower=-500, upper=500)
eval_tv_opt <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=opt$par))

plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_tv_opt, col="red")

#############################################################################
## Possibility contour literature combiners
#############################################################################

eval_hose <- sapply(theta_vec, function(t) hose_pl(t, Y, sigma, n))/
  max(sapply(theta_vec, function(t) hose_pl(t, Y, sigma, n)))
eval_HH <- sapply(theta_vec, function(t) HH_pl(t, Y, sigma, n, weight=1/sigma^2/sum(1/sigma^2) ))/
  max(sapply(theta_vec, function(t) HH_pl(t, Y, sigma, n, weight=1/sigma^2/sum(1/sigma^2) )))


#############################################################################
## p-to-E-to-p calibrator combiners
#############################################################################

eval_e_to_p_AR <- sapply(theta_vec, function(t) e_to_p_AR(t, Y, sigma, n))
eval_e_to_p_kappa <- sapply(theta_vec, function(t) e_to_p_kappa(t, Y, sigma, n, weight=rep(0.52,n)))

#############################################################################
## Make the figures for p-value combiners
#############################################################################

pdf("simple-normal-opt-ind.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
dev.off()

pdf("simple-normal-tippett-fisher.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_tippett, col="red")
lines(theta_vec, eval_fisher, col="green")
dev.off()

pdf("simple-normal-others.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_pearson, col="red")
lines(theta_vec, eval_MG, col="green")
lines(theta_vec, eval_edgington, col="mediumblue")
lines(theta_vec, eval_stouffer, col="purple")
dev.off()

pdf("simple-normal-prod.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_prod, col="red")
lines(theta_vec, eval_prod_transf, col="green")
lines(theta_vec, eval_fisher, col="purple")
dev.off()

pdf("simple-normal-kolmogorov.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
# lines(theta_vec, eval_kolmogorov_min, col="red")
# lines(theta_vec, eval_kolmogorov_max, col="green")
lines(theta_vec, eval_kolmogorov_harmonic, col="red")
lines(theta_vec, eval_kolmogorov_geometric, col="green")
lines(theta_vec, eval_kolmogorov_arithmetic, col="purple")
lines(theta_vec, eval_kolmogorov_neg, col="orange")
lines(theta_vec, eval_kolmogorov_pos, col="pink")
dev.off()

r <- -0.5
eval_kolmogorov_neg <- sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=r))/
  max(sapply(theta_vec, function(t) kolmogorov_pl(t, Y, sigma, n, r=r)))
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_kolmogorov_neg, col="orange")

plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_kolmogorov_BA, col="red")
lines(theta_vec, eval_kolmogorov_BG, col="green")

#############################################################################
## Make the figures for optimization combiners
#############################################################################

pdf("simple-normal-median.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_median, col="red")
lines(theta_vec, eval_p25, col="green")
dev.off()


#############################################################################
## Make the figures for possibility contour literature combiners
#############################################################################

pdf("simple-normal-hose.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_hose_norm, col="red")
# lines(theta_vec, eval_valid_hose, col="green")
dev.off()

pdf("simple-normal-HH.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_hose, col="red")
lines(theta_vec, eval_HH, col="green")
# lines(theta_vec, eval_valid_HH, col="green")
dev.off()

## Both of these look pretty good
v <- 1/c(3,8,2)
pdf("simple-normal-HH_272.pdf", height=5, width=5)
v <- 1/c(2,7,2)
eval_HH <- sapply(theta_vec, function(t) HH_pl(t, Y, sigma, n, weight=v/sum(v) ))
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="blue")
lines(theta_vec, eval_individual2, col="orange")
lines(theta_vec, eval_individual3, col="purple")
lines(theta_vec, eval_HH, col="red")
dev.off()

#############################################################################
## Make the figures for p-to-E-to-p calibrator combiners
#############################################################################

pdf("simple-normal-p_to_e.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, eval_individual1, col="grey")
lines(theta_vec, eval_individual2, col="grey")
lines(theta_vec, eval_individual3, col="grey")
lines(theta_vec, eval_e_to_p_kappa, col="green")
lines(theta_vec, eval_e_to_p_AR, col="red")
dev.off()

