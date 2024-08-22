// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat matrix_inv(const arma::mat& X){
  arma::mat X_inv = inv(X);
  return(X_inv);
}

// [[Rcpp::export]]
arma::vec lognormal_MLE(const arma::vec& Y){
  double mu_hat = mean(log(Y));
  double gamma_hat = mean(pow(log(Y) - mu_hat, 2.0));
  arma::vec theta_MLE = zeros<vec>(2);
  theta_MLE(0) = mu_hat;
  theta_MLE(1) = gamma_hat;
  return(theta_MLE);
}

// [[Rcpp::export]]
arma::mat lognormal_FI(const arma::vec& Y, const double& n, const arma::vec& theta){
  arma::mat J_ni = zeros<mat>(2,2);
  J_ni(0,0) = n/theta(1);
  J_ni(0,1) = sum(log(Y)-theta(0))/pow(theta(1), 2.0);
  J_ni(1,0) = J_ni(0,1);
  J_ni(1,1) = sum(pow(log(Y) - theta(0), 2.0)) / pow(theta(1), 3.0) - n/(2*pow(theta(1), 2.0));
  return(J_ni);
} 

// [[Rcpp::export]]
double pareto_QF(double const& u, double const& theta){
  double QF = pow(1-u, -1/theta);
  return(QF);
}

// [[Rcpp::export]]
double logdlnorm(const arma::vec& Y, const double& n, const arma::vec& theta){
  double ll = -sum(log(Y)) - n*log(theta(1))/2 - sum(pow(log(Y) - theta(0), 2.0)) / (2*theta(1));
  return(ll);
}

// [[Rcpp::export]]
double logdbinom(const arma::vec& Y, const arma::vec& mu){
  double ll = sum(Y % log(mu) + (1-Y)%log(1-mu));
  return(ll);
}

// [[Rcpp::export]]
double lognormal_individual_IM(const arma::vec& Y, const double& n, const double& mu, const double& gamma, const arma::vec& theta_MLE){
  double ll1 = - n*log(gamma)/2 - sum(pow(log(Y) - mu, 2.0)) / (2*gamma);
  double ll2 = - n*log(theta_MLE(1))/2 - sum(pow(log(Y) - theta_MLE(0), 2.0)) / (2*theta_MLE(1));
  double lr = exp(ll1 - ll2);
  return(lr);
}

// [[Rcpp::export]]
double pareto_individual_IM(const arma::vec& Y, const double& n, const double& theta){
  double lr = exp(n*log(theta)+((n/sum(log(Y))+1)-(theta+1))*sum(log(Y)) - n*log(n/sum(log(Y))) );
  return(lr);
}

// [[Rcpp::export]]
double logistic_individual_IM(const arma::vec& mu, const arma::vec& mu_MLE, const arma::vec& Y){
  double lr = exp(logdbinom(Y, mu) - logdbinom(Y, mu_MLE));
  return(lr);
}

// [[Rcpp::export]]
arma::vec valid_lognormal_individual_IM(const arma::mat & theta, const arma::vec& Y, const double& n, const double& M){
  int T = theta.n_rows;
  arma::vec theta_MLE = zeros<vec>(2);
  theta_MLE(0) = mean(log(Y));
  theta_MLE(1) = mean(pow(log(Y) - theta_MLE(0), 2.0));
  arma::vec pl = zeros<vec>(T);
  
  for(int t=0; t<T; t++){
    double observed_contour = lognormal_individual_IM(Y, n, theta(t,0), theta(t,1), theta_MLE);
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::vec Ysim = zeros<vec>(n);
      for(int j=0; j<n; j++){
        Ysim(j) = R::rlnorm(theta(t,0), pow(theta(t,1), 0.5)); 
      }
      arma::vec theta_MLE_sim = zeros<vec>(2);
      theta_MLE_sim(0) = mean(log(Ysim));
      theta_MLE_sim(1) = mean(pow(log(Ysim) - theta_MLE_sim(0), 2.0));
      double sim_contour = lognormal_individual_IM(Ysim, n, theta(t,0), theta(t,1), theta_MLE_sim);
      if(sim_contour <= observed_contour) res(sim) = 1;
    }
    pl(t) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec valid_pareto_individual_IM(const arma::vec & theta, arma::vec const& Y, const double &n, const double& M){
  int T = theta.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  for(int t=0; t<T; t++){
    double observed_contour = pareto_individual_IM(Y, n, theta(t));
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::vec Ysim = zeros<vec>(n);
      for(int j=0; j<n; j++){
        Ysim(j) = pareto_QF(R::runif(0,1), theta(t)); 
      }
      double sim_contour = pareto_individual_IM(Ysim, n, theta(t));
      if(sim_contour <= observed_contour) res(sim) = 1;
    }
    pl(t) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec valid_logistic_individual_IM(const arma::mat & theta, const arma::vec theta_MLE, const arma::vec& Y, const arma::mat covariates,
                                       const double& n, const double& p, const double& M){
  int T = theta.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  arma::mat Xtheta = covariates * theta.t();
  arma::mat mu = exp(Xtheta)/(1+exp(Xtheta));
  
  arma::vec Xtheta_MLE = covariates * theta_MLE;
  arma::vec mu_MLE = exp(Xtheta_MLE)/(1+exp(Xtheta_MLE)); 

  for(int t=0; t<T; t++){
    double observed_contour = logistic_individual_IM(mu.col(t), mu_MLE, Y);
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::vec Ysim = zeros<vec>(n);
      for(int j=0; j<n; j++){
        Ysim(j) = R::rbinom(1, mu(j,t)); 
      }
      
      // GLM optimization modified from package here: github.com/dirkschumacher/rcppglm
      arma::mat Q, R;
      arma::colvec s = arma::zeros<arma::colvec>(p);
      arma::colvec s_old;
      arma::colvec eta = arma::ones<arma::colvec>(n);
      arma::qr_econ(Q, R, covariates);
      
      for (int it = 0; it < 10000; it++) {
        s_old = s;
        const arma::colvec mu = (arma::exp(eta) / (1.0 + arma::exp(eta)));
        const arma::colvec mu_p = arma::exp(eta) / arma::square(arma::exp(eta) + 1.0);
        const arma::colvec z = eta + (Ysim - mu) / mu_p;
        const arma::colvec W = arma::square(mu_p) / (mu % (1.0 - mu));
        const arma::mat C = arma::chol(Q.t() * (Q.each_col() % W));
        const arma::colvec s1 = arma::solve(arma::trimatl(C.t()), Q.t() * (W % z));
        s = arma::solve(arma::trimatu(C), s1);
        eta = Q * s;

        const bool is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < 0.000001;
        if (is_converged) break;
      }
      
      arma::vec theta_MLE_sim = arma::solve(arma::trimatu(R), Q.t() * eta);
      arma::vec Xtheta_MLE_sim = covariates * theta_MLE_sim;
      arma::vec mu_MLE_sim = exp(Xtheta_MLE_sim)/(1+exp(Xtheta_MLE_sim));
      double sim_contour = logistic_individual_IM(mu.col(t), mu_MLE_sim, Ysim);
      if(sim_contour <= observed_contour) res(sim) = 1;
    }
    pl(t) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
double lognormal_work_likeli(arma::vec const& theta, arma::mat const& theta_MLE, arma::field<arma::mat> const& J_ni, int const& k){
  arma::mat s_J = zeros<mat>(2,2);
  arma::vec main = zeros<vec>(2);
  for(int i=0; i<k; i++){
    s_J = s_J + J_ni(i);
    main += J_ni(i) * theta_MLE.row(i).t();
  }
  arma::mat s_J_inv = inv(s_J);
  arma::vec theta_n = s_J_inv * main;
  double ll = exp(arma::as_scalar( - (theta_n-theta).t() * s_J * (theta_n - theta)) / 2);
  return(ll);
}

// [[Rcpp::export]]
double pareto_work_likeli(double const& theta, const arma::vec& theta_MLE, arma::vec const& J_ni, int const& k){
  double s_J = sum(J_ni);
  double theta_n = sum(theta_MLE % J_ni) / s_J;
  double ll = exp(-pow(theta_n - theta, 2.0) * s_J / 2);
  return(ll);
}

// [[Rcpp::export]]
double logistic_work_likeli(arma::vec const& theta, arma::mat const& theta_MLE, arma::field<arma::mat> const& J_ni, int const& p, int const& k){
  arma::mat s_J = zeros<mat>(p,p);
  arma::vec main = zeros<vec>(p);
  for(int i=0; i<k; i++){
    s_J = s_J + J_ni(i);
    main += J_ni(i) * theta_MLE.row(i).t();
  }
  arma::mat s_J_inv = inv(s_J);
  arma::vec theta_n = s_J_inv * main;
  double ll = exp(arma::as_scalar( - (theta_n-theta).t() * s_J * (theta_n - theta)) / 2);
  return(ll);
}


// [[Rcpp::export]]
arma::vec lognormal_pl_mid_exact(arma::mat const& theta, arma::mat const& theta_MLE, int const& k, 
                       arma::field<arma::mat> const& J_ni, arma::vec const& n_i, int const& M){
  double T = theta.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  for(int l=0; l<T; l++){
    double wl_obs = lognormal_work_likeli(theta.row(l).t(), theta_MLE, J_ni, k);
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,2);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = R::rlnorm(theta(l,0), pow(theta(l,1), 0.5)); 
        }
        theta_MLE_sim(i,0) = mean(log(Ysim_i));
        theta_MLE_sim(i,1) = mean(pow(log(Ysim_i) - theta_MLE_sim(i,0), 2.0));
        J_ni_sim(i) = lognormal_FI(Ysim_i, n_i(i), theta_MLE_sim.row(i).t()); 
      }
      double wl_sim = lognormal_work_likeli(theta.row(l).t(), theta_MLE_sim, J_ni_sim, k);
      if(wl_sim <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec pareto_pl_mid_exact(arma::vec const& theta, arma::vec const& theta_MLE, int const& k, 
                              arma::vec const& J_ni, arma::vec const& n_i, int const& M){
  double T = theta.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  for(int l=0; l<T; l++){
    double wl_obs = pareto_work_likeli(theta(l), theta_MLE, J_ni, k);
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::vec theta_MLE_sim = zeros<vec>(k);
      arma::vec J_ni_sim = zeros<vec>(k);
      for(int i=0; i<k; i++){
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = pareto_QF(R::runif(0,1), theta(l)); 
        }
        theta_MLE_sim(i) = n_i(i)/sum(log(Ysim_i));
        J_ni_sim(i) = n_i(i) / pow(theta_MLE_sim(i),2.0); 
      }
      double wl_sim = pareto_work_likeli(theta(l), theta_MLE_sim, J_ni_sim, k);
      if(wl_sim <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec logistic_pl_mid_exact(arma::mat const& theta, arma::mat const& theta_MLE, int const& k, 
                                arma::field<arma::mat> const& J_ni, arma::vec const& n_i, 
                                arma::field<arma::mat> const& covariates, int const& M){
  int p = theta.n_cols;
  double T = theta.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  arma::field<arma::mat> mu(k);
  for(int i=0; i<k; i++){
    arma::mat Xtheta = covariates(i) * theta.t();
    mu(i) = exp(Xtheta)/(1+exp(Xtheta)); 
  }
  
  for(int l=0; l<T; l++){
    double wl_obs = logistic_work_likeli(theta.row(l).t(), theta_MLE, J_ni, p, k);
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,2);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = R::rbinom(1, mu(i)(j,l)); 
        }
        arma::vec theta_MLE_sim_i;
        if( sum(Ysim_i==0) >= n_i(i)-1 || sum(Ysim_i==1) >= n_i(i)-1 ) {
          theta_MLE_sim_i = zeros<vec>(p);
          theta_MLE_sim.row(i) = theta_MLE_sim_i.t();
          J_ni_sim(i) = zeros<mat>(p,p);
        } else {
          // GLM optimization modified from package here: github.com/dirkschumacher/rcppglm
          arma::mat Q, R, C, C2;
          arma::colvec s = arma::zeros<arma::colvec>(p);
          arma::colvec s1, s_old, mu, mu_p, z, W;
          arma::colvec eta = covariates(i) * theta.row(l).t();
          arma::qr_econ(Q, R, covariates(i));
          bool is_converged, error;
          error = false;
          
          for (int it = 0; it < 10000; it++) {
            s_old = s;
            mu = (arma::exp(eta) / (1.0 + arma::exp(eta)));
            mu_p = arma::exp(eta) / arma::square(arma::exp(eta) + 1.0);
            z = eta + (Ysim_i - mu) / mu_p;
            W = arma::square(mu_p) / (mu % (1.0 - mu));
            C2 = Q.t() * (Q.each_col() % W);
            if(C2.is_symmetric()){
              C = arma::chol(C2);
              s1 = arma::solve(arma::trimatl(C.t()), Q.t() * (W % z));
              s = arma::solve(arma::trimatu(C), s1);
              eta = Q * s;
              is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < 0.000001;
            } else{
              eta = zeros<vec>(n_i(i));
              error = true;
            }
            if (is_converged || error) break;
          }
          
          theta_MLE_sim_i = arma::solve(arma::trimatu(R), Q.t() * eta);
          theta_MLE_sim.row(i) = theta_MLE_sim_i.t();
          J_ni_sim(i) = zeros<mat>(p,p);
          for(int j=0; j<n_i(i); j++){
            double mu_ij = arma::as_scalar(exp(covariates(i).row(j)*theta_MLE_sim_i));
            J_ni_sim(i) = J_ni_sim(i) + covariates(i).row(j).t()*covariates(i).row(j) * mu_ij / pow(1+mu_ij, 2.0);
          } 
        }
      }
      double wl_sim = logistic_work_likeli(theta.row(l).t(), theta_MLE_sim, J_ni_sim, p, k);
      if(wl_sim <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec logistic_pl_mid(arma::mat const& theta, arma::mat const& theta_dag, arma::mat const& theta_MLE, int const& k, 
                                arma::field<arma::mat> const& J_ni, arma::vec const& n_i, 
                                arma::field<arma::mat> const& covariates, int const& M){
  int p = theta.n_cols;
  double T = theta.n_rows;
  double D = theta_dag.n_rows;
  arma::vec pl = zeros<vec>(T);
  
  arma::field<arma::mat> mu_dag(k);
  arma::field<arma::mat> Ysim(k);
  for(int i=0; i<k; i++){
    arma::mat Xtheta_dag = covariates(i) * theta_dag.t();
    mu_dag(i) = exp(Xtheta_dag)/(1+exp(Xtheta_dag)); 
    
    arma::mat Ysim_i = zeros<mat>(n_i(i)*M,D);
    for(int d=0; d<D; d++){
      for(int sim=0; sim<M; sim++){
        for(int j=0; j<n_i(i); j++){
          Ysim_i(sim*n_i(i)+j,d) = R::rbinom(1, mu_dag(i)(j,d));  
        }
      }
    }
    Ysim(i) = Ysim_i;
  }
  
  arma::mat distance = zeros<mat>(T,D);
  for(int d=0; d<D; d++){
    for(int t=0; t<T; t++){
      distance(t,d) = pow(arma::as_scalar( theta.row(t) * theta_dag.row(d).t() ), 0.5);
    }
  }
  
  for(int l=0; l<T; l++){
    double wl_obs = logistic_work_likeli(theta.row(l).t(), theta_MLE, J_ni, p, k);
    
    int closest = distance.row(l).index_min();
    arma::field<arma::vec> mu(k);
    for(int i=0; i<k; i++){
      arma::vec Xtheta = covariates(i) * theta.row(l).t();
      mu(i) = exp(Xtheta)/(1+exp(Xtheta)); 
    }
    
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,2);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::vec Ysim_i = Ysim(i)(span(sim*n_i(i), (sim+1)*n_i(i)-1), closest);
        
        // GLM optimization modified from package here: github.com/dirkschumacher/rcppglm
        arma::mat Q, R, C, C2;
        arma::colvec s = arma::zeros<arma::colvec>(p);
        arma::colvec s1, s_old, mu, mu_p, z, W;
        arma::colvec eta = covariates(i) * theta.row(l).t();
        arma::qr_econ(Q, R, covariates(i));
        bool is_converged, error;
        error = false;
        
        for (int it = 0; it < 10000; it++) {
          s_old = s;
          mu = (arma::exp(eta) / (1.0 + arma::exp(eta)));
          mu_p = arma::exp(eta) / arma::square(arma::exp(eta) + 1.0);
          z = eta + (Ysim_i - mu) / mu_p;
          W = arma::square(mu_p) / (mu % (1.0 - mu));
          C2 = Q.t() * (Q.each_col() % W);
          if(C2.is_symmetric(1e-3)){
            C = arma::chol(C2);
            s1 = arma::solve(arma::trimatl(C.t()), Q.t() * (W % z));
            s = arma::solve(arma::trimatu(C), s1);
            eta = Q * s;
            is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < 0.000001;
          } else{
            eta = zeros<vec>(n_i(i));
            error = true;
          }
          if (is_converged || error) break;
        }
        
        arma::vec theta_MLE_sim_i = arma::solve(arma::trimatu(R), Q.t() * eta);
        theta_MLE_sim.row(i) = theta_MLE_sim_i.t();
        J_ni_sim(i) = zeros<mat>(p,p);
        for(int j=0; j<n_i(i); j++){
          double mu_ij = arma::as_scalar(exp(covariates(i).row(j)*theta_MLE_sim_i));
          J_ni_sim(i) = J_ni_sim(i) + covariates(i).row(j).t()*covariates(i).row(j) * mu_ij / pow(1+mu_ij, 2.0);
        }
      }
      double wl_sim = logistic_work_likeli(theta.row(l).t(), theta_MLE_sim, J_ni_sim, p, k);
      if(wl_sim <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

