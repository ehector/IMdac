# ALL <- paste0("bernoulli-onedag-sim",1:504,".RData")
# missing <- c()
# for(a in 1:504){
#   if(!(ALL[a] %in% Sys.glob("bernoulli-onedag-sim*"))){
#     missing <- c(missing,a)
#   }
# }
# 
# IM_matrix <- naive_IM_matrix <- exact_IM_matrix <- matrix(0, 160001, 504)
# theta_vec_matrix <- list()
# for(sim in setdiff(1:504,missing)){
#   load(paste0("bernoulli-onedag-sim",sim,".RData"))
#   theta_vec_matrix[[sim]] <- theta_vec
#   IM_matrix[,sim] <- middle_IM
#   naive_IM_matrix[,sim] <- naive_IM
# }
# 
# nsim <- 504-length(missing)
# bad_sims <- missing
# IM_matrix <- IM_matrix[,-bad_sims][,1:nsim]
# naive_IM_matrix <- naive_IM_matrix[,-bad_sims][,1:nsim]
# theta_vec_matrix <- theta_vec_matrix[-bad_sims][1:nsim]
# save.image("bernoulli-onedag-inference.RData")

load("/Users/ehector/Dropbox/Projects/IMfusion/simulations/RData files/bernoulli-onedag-inference.RData")

alpha_vec <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
IM_coverages <- naive_coverages <- rep(0, length(alpha_vec))
for(a in 1:length(alpha_vec)){
  alpha <- alpha_vec[a]
  IM_coverages[a] <-
    mean(sapply(1:nsim, function(sim){
      theta_CI <- theta_vec_matrix[[sim]][which(IM_matrix[,sim] > alpha),]
      # lb <- c(min(theta_CI[,1]), min(theta_CI[,2]))
      # ub <- c(max(theta_CI[,1]), max(theta_CI[,2]))
      # as.numeric(all(theta >= lb) && all(theta <= ub))
      as.numeric( sum(theta_CI[1,] == theta)==2 )
    }))
  
  naive_coverages[a] <-
    mean(sapply(1:nsim, function(sim){
      theta_CI <- theta_vec_matrix[[sim]][which(naive_IM_matrix[,sim] > alpha),]
      # lb <- c(min(theta_CI[,1]), min(theta_CI[,2]))
      # ub <- c(max(theta_CI[,1]), max(theta_CI[,2]))
      as.numeric( sum(theta_CI[1,] == theta)==2 )
    }))
}
rbind(rev(1-alpha_vec), rev(IM_coverages), rev(naive_coverages))

paste(rev(100*(1-alpha_vec)), collapse=" & ")
paste(rev(naive_coverages)*100, collapse=" & ")
paste(rev(IM_coverages)*100, collapse=" & ")

## Calculate confidence interval lengths
alpha_vec <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
IM_lengths <- naive_lengths <- rep(0, length(alpha_vec))
for(a in 1:length(alpha_vec)){
  alpha <- alpha_vec[a]
  IM_lengths[a] <-
    mean(sapply(1:nsim, function(sim){
      theta_CI <- theta_vec_matrix[[sim]][which(IM_matrix[,sim] > alpha),]
      # lb <- c(min(theta_CI[,1]), min(theta_CI[,2]))
      # ub <- c(max(theta_CI[,1]), max(theta_CI[,2]))
      # ub - lb
      mm <- which.max(IM_matrix[,sim])
      max(abs(theta_CI[,1] - theta_vec_matrix[[sim]][mm,1])) * max(abs(theta_CI[,2] - theta_vec_matrix[[sim]][mm,2])) * pi
    }))
  naive_lengths[a] <-
    mean(sapply(1:nsim, function(sim){
      theta_CI <- theta_vec_matrix[[sim]][which(naive_IM_matrix[,sim] > alpha),]
      # lb <- c(min(theta_CI[,1]), min(theta_CI[,2]))
      # ub <- c(max(theta_CI[,1]), max(theta_CI[,2]))
      # ub-lb
      mm <- which.max(naive_IM_matrix[,sim])
      max(abs(theta_CI[,1] - theta_vec_matrix[[sim]][mm,1])) * max(abs(theta_CI[,2] - theta_vec_matrix[[sim]][mm,2])) * pi
    }))
}
paste(rev(100*(1-alpha_vec)), collapse=" & ")
paste(rev(signif(naive_lengths*10,3)), collapse=" & ")
paste(rev(signif(IM_lengths*10,3)), collapse=" & ")

