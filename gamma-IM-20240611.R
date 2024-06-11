set.seed(600)
theta <- 0.5
k <- 3
n_i <- c(5,10,15)
Y <- sapply(1:k, function(i) rgamma(1, n_i[i], rate=theta) )

expo_individual_IM <- function(theta, Y){
  z <- theta*Y*exp(-theta*Y)
  1-exp(lamW::lambertW0(-z)) + exp(lamW::lambertWm1(-z))
}
gamma_individual_IM <- function(theta, Y, n){
  z <- (theta*Y)^n * exp(-theta*Y)
  expint::gammainc(n, -n*lamW::lambertWm1(-z^(1/n)/n))/factorial(n-1) -
    expint::gammainc(n, -n*lamW::lambertW0(-z^(1/n)/n))/factorial(n-1) + 1
}
gamma_combined_IM <- function(theta, Y, n){
  z <- (theta*sum(Y))^n * exp(-theta*sum(Y))
  expint::gammainc(n, -n*lamW::lambertWm1(-z^(1/n)/n))/factorial(n-1) -
    expint::gammainc(n, -n*lamW::lambertW0(-z^(1/n)/n))/factorial(n-1) + 1
}
recreate_gamma_combined_IM <- function(t, Y, n_i, k, individual_IMs, MC_samples_sort){
  G_inv_pi <-
    sapply(1:k, function(i){
      pi <- individual_IMs[[i]][t]
      if(pi<1e-4) return(0)
      else return(MC_samples_sort[[i]][num_samp*pi])
    })
  theta_Y <- sapply(1:k, function(i) {
    if(Y[i] >= n_i[i]/theta_vec[t]){
      return( -n_i[i] * lamW::lambertWm1(-G_inv_pi[i]^(1/n_i[i]) / n_i[i]) )
    }
    else return( -n_i[i] * lamW::lambertW0(-G_inv_pi[i]^(1/n_i[i]) / n_i[i]) )
  })
  z <- sum(theta_Y)^sum(n_i) * exp(-sum(theta_Y))
  prob <- 1 +
    expint::gammainc(sum(n_i), -sum(n_i)*lamW::lambertWm1( -z^(1/sum(n_i))/sum(n_i)) )/factorial(sum(n_i)-1) - 
    expint::gammainc(sum(n_i), -sum(n_i)*lamW::lambertW0(  -z^(1/sum(n_i))/sum(n_i)) )/factorial(sum(n_i)-1)
  prob
}

theta_vec <- seq(1e-10,1.5,length=1000)

individual_IMs <- lapply(1:k, function(i) sapply(theta_vec, function(t) gamma_individual_IM(theta=t, Y=Y[i], n=n_i[i]) ))
eval_optimal <- sapply(theta_vec, function(t) gamma_combined_IM(theta=t, Y=Y, n=sum(n_i)) )

pdf("simple-gamma-opt-ind.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)))
lines(theta_vec, individual_IMs[[1]], col="grey")
lines(theta_vec, individual_IMs[[2]], col="grey")
lines(theta_vec, individual_IMs[[3]], col="grey")
dev.off()

num_samp <- 10000
MC_samples_ind <- lapply(1:k, function(i) rgamma(num_samp, shape=n_i[i], rate=1) )
MC_samples_sort <- lapply(1:k, function(i) sort(MC_samples_ind[[i]]^(n_i[i]) * exp(-MC_samples_ind[[i]])))
eval_optimal_reconstruct <- sapply(1:length(theta_vec), function(t) recreate_gamma_combined_IM(t, Y, n_i, k, individual_IMs, MC_samples_sort) )
plot(theta_vec, eval_optimal_reconstruct, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)), ylim=c(0,1))
lines(theta_vec, eval_optimal, col="gray")

eval_optimal_approx <- sapply(1:length(theta_vec), function(t) {
  sigma <- theta_vec[t]/sqrt(n_i)
  1-abs(2*pnorm( (sum(1/sigma^2))^(-1/2) * sum(-sign(Y-theta_vec[t])/sigma * qnorm(sapply(individual_IMs, function(x) x[t])/2)) ) - 1)
})
eval_fisher <- sapply(1:length(theta_vec), function(t) {
  1-pchisq(-2*sum(log(sapply(individual_IMs, function(x) x[[t]]))), df=2*k)
})

pdf("simple-gamma-opt-ind-approx.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)), ylim=c(0,1))
lines(theta_vec, eval_optimal_approx/max(eval_optimal_approx), col="red")
lines(theta_vec, eval_fisher/max(eval_fisher), col="green")
lines(theta_vec, individual_IMs[[1]], col="grey")
lines(theta_vec, individual_IMs[[2]], col="grey")
lines(theta_vec, individual_IMs[[3]], col="grey")
dev.off()


theta_hat <- n_i/Y
variances <- theta_hat^2/n_i
theta_tilde <- sum(theta_hat / variances)/sum(1/variances)
eval_optimal_approx2 <- sapply(theta_vec, function(t) optimal_pl(t, theta_tilde, sqrt(variances), n=k))

pdf("simple-gamma-opt-ind-approx2.pdf", height=5, width=5)
plot(theta_vec, eval_optimal, type="l", xlab=expression(theta), ylab=expression(pi[y](theta)), ylim=c(0,1))
lines(theta_vec, eval_optimal_approx2, col="red")
lines(theta_vec, eval_fisher/max(eval_fisher), col="green")
lines(theta_vec, individual_IMs[[1]], col="grey")
lines(theta_vec, individual_IMs[[2]], col="grey")
lines(theta_vec, individual_IMs[[3]], col="grey")
dev.off()
