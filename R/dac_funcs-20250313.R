create_rng <- function(seed = NULL) {
  # Set the random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create a random number generator object
  rng <- new.env(parent = emptyenv())
  
  rng$g_and_k <- function(size, A, B, g, k, c=0.8) {
    z <- rnorm(size, mean = 0, sd = 1)
    zeros <- rep(0, size)
    z <- z + zeros
    A <- A + zeros
    B <- B + zeros
    g <- g + zeros
    k <- k + zeros
    c <- c + zeros
    z_squared <- z^2
    term1 <- (1 + c * tanh(g * z/2))
    term2 <- z * (1 + z_squared)^k
    term1[g == 0] <- 1
    zbig <- which(is.infinite(z_squared))
    term2[zbig] <- sign(z[zbig]) * abs(z[zbig])^(1 + 2 * k[zbig])
    sample <- matrix(A + B * term1 * term2, ncol=1)
    return(sample)
  }
  
  rng$g_and_k_PM25 <- function(size, A, sigma_params, g, k, c=0.8, covariates_sigma) {
    z <- rnorm(size, mean = 0, sd = 1)
    zeros <- rep(0, size)
    z <- z + zeros
    A <- A + zeros
    B <- exp(covariates_sigma %*% sigma_params)
    g <- g + zeros
    k <- k + zeros
    c <- c + zeros
    z_squared <- z^2
    term1 <- (1 + c * tanh(g * z/2))
    term2 <- z * (1 + z_squared)^k
    term1[g == 0] <- 1
    zbig <- which(is.infinite(z_squared))
    term2[zbig] <- sign(z[zbig]) * abs(z[zbig])^(1 + 2 * k[zbig])
    sample <- matrix(A + B * term1 * term2, ncol=1)
    return(sample)
  }
  
  rng$alpha_stable <- function(size, alpha, beta, c, mu){
    U <- runif(size, -pi/2, pi/2)
    W <- rexp(size, 1)
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
    return(matrix(Y, ncol=1))
  }
  
  rng$uniform <- function(size, min = 0, max = 1) {
    runif(size, min=min, max=max)
  }
  
  return(rng)
}

train_g_and_k_emulator <- function(n, num_summary_dim){
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-20, max=20), RNG$uniform(size = 1L, min=0, max=20),
             RNG$uniform(size = 1L, min=-5, max=5), RNG$uniform(size = 1L, min=-1/2, max=5)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("A","B","g","k"))
  
  likelihood_fun <- function(params, n_obs = n) {
    return(RNG$g_and_k(size = n_obs, A=params[1], B=params[2], g=params[3], k=params[4]))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = 4L,
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  num_epochs <- 1000L
  num_its <- 1000L
  history <- trainer$train_online(
    epochs = num_epochs,
    iterations_per_epoch = num_its,
    batch_size = 32L,
    validation_sims = 200L,
    use_act_norm = TRUE
  )
  
  return(amortizer)
}

load_trained_g_and_k_emulator <- function(n, num_summary_dim, path){
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-20, max=20), RNG$uniform(size = 1L, min=0, max=20),
             RNG$uniform(size = 1L, min=-5, max=5), RNG$uniform(size = 1L, min=-1/2, max=5)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("A","B","g","k"))
  
  likelihood_fun <- function(params, n_obs = n) {
    return(RNG$g_and_k(size = n_obs, A=params[1], B=params[2], g=params[3], k=params[4]))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = 4L,
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  amortizer$load_weights(path)
  
  return(amortizer)
}

train_g_and_k_emulator_PM25 <- function(n, num_summary_dim, covariates_sigma){
  p_sigma <- ncol(covariates_sigma)
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-20, max=20), 
             sapply(1:p_sigma, function(q) RNG$uniform(size = 1L, min=-2, max=2)),
             RNG$uniform(size = 1L, min=-5, max=5),
             RNG$uniform(size = 1L, min=-1/2, max=5)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("A", paste0("sigma_",1:p_sigma),"g","k"))
  
  likelihood_fun <- function(params, n_obs = n) {
    simulations <- RNG$g_and_k_PM25(size = n_obs, A=params[1], sigma_params=params[2:(p_sigma+1)], 
                                    g=params[p_sigma+2], k=params[p_sigma+3], 
                                    covariates_sigma=covariates_sigma)
    return(cbind(simulations, covariates_sigma))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = as.integer(p_sigma+3),
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  num_epochs <- 5000L
  num_its <- 1000L
  history <- trainer$train_online(
    epochs = num_epochs,
    iterations_per_epoch = num_its,
    batch_size = 32L,
    validation_sims = 200L,
    use_act_norm = TRUE
  )
  
  return(amortizer)
}

load_trained_g_and_k_emulator_PM25 <- function(n, num_summary_dim, covariates_sigma, path){
  p_sigma <- ncol(covariates_sigma)
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-20, max=20), 
             sapply(1:p_sigma, function(q) RNG$uniform(size = 1L, min=-2, max=2)),
             RNG$uniform(size = 1L, min=-5, max=5),
             RNG$uniform(size = 1L, min=-1/2, max=5)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("A", paste0("sigma_",1:p_sigma),"g","k"))
  
  likelihood_fun <- function(params, n_obs = n) {
    simulations <- RNG$g_and_k_PM25(size = n_obs, A=params[1], sigma_params=params[2:(p_sigma+1)], 
                                    g=params[p_sigma+2], k=params[p_sigma+3], 
                                    covariates_sigma=covariates_sigma)
    return(cbind(simulations, covariates_sigma))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = as.integer(p_sigma+3),
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  amortizer$load_weights(path)
  
  return(amortizer)
}

train_alpha_stable_emulator <- function(n, num_summary_dim, alpha){
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-1, max=1),
             RNG$uniform(size = 1L, min=0, max=10), 
             RNG$uniform(size = 1L, min=-20, max=20)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("beta","c","mu"))
  
  likelihood_fun <- function(params, n_obs = n) {
    return(RNG$alpha_stable(size = n_obs, alpha=alpha, beta=params[1], c=params[2], mu=params[3]))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = 3L,
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  num_epochs <- 1000L
  num_its <- 1000L
  history <- trainer$train_online(
    epochs = num_epochs,
    iterations_per_epoch = num_its,
    batch_size = 32L,
    validation_sims = 200L,
    use_act_norm = TRUE
  )
  
  return(amortizer)
}

load_trained_alpha_stable_emulator <- function(n, num_summary_dim, alpha, path){
  RNG <- create_rng(seed = 2023)
  
  prior_fun <- function() {
    return(c(RNG$uniform(size = 1L, min=-1, max=1),
             RNG$uniform(size = 1L, min=0, max=10), 
             RNG$uniform(size = 1L, min=-20, max=20)))
  }
  
  # Define Prior using BayesFlow's Prior class
  prior <- bf$simulation$Prior(prior_fun = prior_fun, param_names=c("beta","c","mu"))
  
  likelihood_fun <- function(params, n_obs = n) {
    return(RNG$alpha_stable(size = n_obs, alpha=alpha, beta=params[1], c=params[2], mu=params[3]))
  }
  
  # Define Simulator using BayesFlow's Simulator class
  simulator <- bf$simulation$Simulator(simulator_fun = likelihood_fun)
  
  model <- bf$simulation$GenerativeModel(prior=prior, simulator=simulator)
  
  summary_net <- bf$networks$DeepSet(summary_dim = num_summary_dim)
  
  inference_net <- bf$networks$InvertibleNetwork(
    num_params = 3L,
    num_coupling_layers = 6L,
    coupling_settings = list(
      "dense_args" = list(kernel_regularizer = NULL),
      "dropout" = FALSE
    ),
    permutation="learnable"
  )
  
  amortizer <- bf$amortizers$AmortizedPosterior(inference_net, summary_net)
  
  trainer <- bf$trainers$Trainer(amortizer = amortizer, generative_model = model)
  
  amortizer$load_weights(path)
  
  return(amortizer)
}

g_and_k_pl_mid_all_profile_NN <- function(theta, theta_dag, large_n_estimate, k, s_J, n_i, p, M=1000, n_samples, amortizer){
  T_ <- nrow(theta)
  pl <- matrix(NA,T_,p)
  wl_sim <- matrix(NA,M,p)
  
  Ysim <- lapply(1:k, function(i) {
    yy <- matrix(rgk(M*n_i[i], A = theta_dag[1], B = theta_dag[2], g = theta_dag[3], k = theta_dag[4], c=0.8), nrow=M)
    return(yy)
  })
  
  for(m in 1:M){
    post_samples <- lapply(1:k, function(i){
      conf_data <- list()
      conf_data$summary_conditions <- array(Ysim[[i]][m,], dim=c(1,n_i[i],1))
      conf_data$direct_conditions <- NULL
      post_samples <- amortizer$sample(conf_data, n_samples=n_samples)
      MLE <- colMeans(post_samples)
      var <- var(post_samples)
      return(list(MLE=MLE, var=var))
    })
    theta_MLE_sim <- t(sapply(1:k, function(i) post_samples[[i]]$MLE))
    J_ni_sim <- lapply(1:k, function(i) solve(post_samples[[i]]$var))
    
    s_J_sim <- Reduce("+", J_ni_sim)
    large_n_estimate_sim <- drop(solve(s_J_sim) %*% Reduce("+", lapply(1:k, function(i) J_ni_sim[[i]]%*%theta_MLE_sim[i,]) ))
    
    wl_sim[m,] <- sapply(1:p, function(q) exp(drop( - (large_n_estimate_sim[q] - theta_dag[q])^2 * s_J_sim[q,q]) / 2) )
  }
  
  for(l in 1:T_){
    pl[l,] <- sapply(1:p, function(q) {
      wl_obs <- exp(drop( - (large_n_estimate[q] - theta[l,q])^2 * s_J[q,q]) / 2)
      mean(wl_sim[,q] <= wl_obs, na.rm=T)
      })
  }
  
  return(pl)
}

g_and_k_pl_mid_all_profile_NN_PM25 <- function(theta, theta_dag, large_n_estimate, k, s_J, n_i, p, M=1000, n_samples, amortizer_list, covariates_sigma){
  T_ <- nrow(theta)
  pl <- matrix(NA,T_,p)
  wl_sim <- matrix(NA,M,p)
  
  Ysim <- lapply(1:k, function(i) {
    yy <- matrix(rgk(M*n_i[i], A = theta_dag[1], B = exp(drop(covariates_sigma[[i]]%*%theta_dag[2:9])), 
                     g = theta_dag[10], k = theta_dag[11], c=0.8), nrow=M)
    return(yy)
  })
  
  for(m in 1:M){
    post_samples <- lapply(1:k, function(i){
      conf_data <- list()
      conf_data$summary_conditions <- array(cbind(Ysim[[i]][m,], covariates_sigma[[i]]), dim=c(1,n_i[i],1+ncol(covariates_sigma[[i]])) )
      conf_data$direct_conditions <- NULL
      post_samples <- amortizer_list[[i]]$sample(conf_data, n_samples=n_samples)
      MLE <- colMeans(post_samples)
      var <- var(post_samples)
      return(list(MLE=MLE, var=var))
    })
    theta_MLE_sim <- t(sapply(1:k, function(i) post_samples[[i]]$MLE))
    J_ni_sim <- lapply(1:k, function(i) solve(post_samples[[i]]$var))
    
    s_J_sim <- Reduce("+", J_ni_sim)
    large_n_estimate_sim <- drop(solve(s_J_sim) %*% Reduce("+", lapply(1:k, function(i) J_ni_sim[[i]]%*%theta_MLE_sim[i,]) ))
    
    wl_sim[m,] <- sapply(1:p, function(q) exp(drop( - (large_n_estimate_sim[q] - theta_dag[q])^2 * s_J_sim[q,q]) / 2) )
  }
  
  for(l in 1:T_){
    pl[l,] <- sapply(1:p, function(q) {
      wl_obs <- exp(drop( - (large_n_estimate[q] - theta[l,q])^2 * s_J[q,q]) / 2)
      mean(wl_sim[,q] <= wl_obs, na.rm=T)
    })
  }
  
  return(pl)
}

alpha_stable_pl_mid_all_profile_NN <- function(theta, alpha, theta_dag, large_n_estimate, k, s_J, n_i, p, M=1000, n_samples, amortizer){
  T_ <- nrow(theta)
  pl <- matrix(NA,T_,p)
  wl_sim <- matrix(NA,M,p)
  
  Ysim <- lapply(1:k, function(i) {
    U <- runif(M*n_i[i], -pi/2, pi/2)
    W <- rexp(M*n_i[i], 1)
    zeta <- -theta_dag[1] * tan(pi*alpha/2)
    if(alpha==1){
      xi <- pi/2
      X <- 1/xi*( ( pi/2 + theta_dag[1]*U)*tan(U) - theta_dag[1]*log( (pi/2*W*cos(U))/(pi/2+theta_dag[1]*U) ))
      Y <- theta_dag[2]*X + 2*theta_dag[1]*theta_dag[2]*log(theta_dag[2])/pi + theta_dag[3]
    } else{
      xi <- atan(-zeta)/alpha
      X <- (1+zeta^2)^(1/(2*alpha)) * 
        (sin(alpha*(U+zeta)))/(cos(U))^(1/alpha) * 
        ( (cos(U - alpha*(U+xi)))/W)^((1-alpha)/alpha)
      Y <- theta_dag[2]*X + theta_dag[3]
    }
    yy <- matrix(Y, nrow=M)
    return(yy)
  })
  
  for(m in 1:M){
    post_samples <- lapply(1:k, function(i){
      conf_data <- list()
      conf_data$summary_conditions <- array(Ysim[[i]][m,], dim=c(1,n_i[i],1))
      conf_data$direct_conditions <- NULL
      post_samples <- amortizer$sample(conf_data, n_samples=n_samples)
      MLE <- colMeans(post_samples)
      var <- var(post_samples)
      return(list(MLE=MLE, var=var))
    })
    theta_MLE_sim <- t(sapply(1:k, function(i) post_samples[[i]]$MLE))
    J_ni_sim <- lapply(1:k, function(i) solve(post_samples[[i]]$var))
    
    s_J_sim <- Reduce("+", J_ni_sim)
    large_n_estimate_sim <- drop(solve(s_J_sim) %*% Reduce("+", lapply(1:k, function(i) J_ni_sim[[i]]%*%theta_MLE_sim[i,]) ))
    
    wl_sim[m,] <- sapply(1:p, function(q) exp(drop( - (large_n_estimate_sim[q] - theta_dag[q])^2 * s_J_sim[q,q]) / 2) )
  }
  
  for(l in 1:T_){
    pl[l,] <- sapply(1:p, function(q) {
      wl_obs <- exp(drop( - (large_n_estimate[q] - theta[l,q])^2 * s_J[q,q]) / 2)
      mean(wl_sim[,q] <= wl_obs, na.rm=T)
    })
  }
  
  return(pl)
}
