// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
double logdpois(const arma::vec& Y, const arma::vec& mu){
  double ll = sum(Y % log(mu) - mu);
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
double poisson_individual_IM(const arma::vec& mu, const arma::vec& mu_MLE, const arma::vec& Y){
  double lr = exp(logdpois(Y, mu) - logdpois(Y, mu_MLE));
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
double exponential_work_likeli(double const& theta, arma::vec const& theta_MLE, arma::vec const& J_ni, int const& k){
  double s_J = sum(J_ni);
  double theta_n = arma::as_scalar(J_ni.t() * theta_MLE)/s_J;
  double ll = exp(arma::as_scalar( - s_J * pow(theta_n - theta, 2.0)) / 2);
  return(ll);
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
double poisson_work_likeli(arma::vec const& theta, arma::mat const& theta_MLE, arma::field<arma::mat> const& J_ni, int const& p, int const& k){
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
double poisson_profile_work_likeli(double const& theta, double const& meta_estimate, double const& s_J){
  double ll = exp(arma::as_scalar( - pow(meta_estimate - theta, 2.0) * s_J) / 2);
  return(ll);
}

// [[Rcpp::export]]
arma::vec exponential_pl_mid_P(arma::vec const& theta, arma::vec const& theta_dag, arma::vec const& theta_MLE, int const& k, 
                             arma::vec const& J_ni, arma::vec const& n_i, int const& M){
  double T = theta.size();
  double D = theta_dag.size();
  arma::vec distance = zeros<vec>(D);
  arma::mat wl_sim = zeros<mat>(D,M);
  arma::vec pl = zeros<vec>(T);
  
  for(int d=0; d<D; d++){
    for(int sim=0; sim<M; sim++){
      arma::vec theta_MLE_sim = zeros<vec>(k);
      arma::vec J_ni_sim = zeros<vec>(k);
      for(int i=0; i<k; i++){
        double Ysim_i = R::rgamma(n_i(i), 1/theta_dag(d));
        
        theta_MLE_sim(i) = n_i(i)/Ysim_i;
        J_ni_sim(i) = n_i(i)/pow(theta_MLE_sim(i), 2.0); 
      }
      wl_sim(d,sim) = exponential_work_likeli(theta_dag(d), theta_MLE_sim, J_ni_sim, k);
    }
  }
  
  for(int l=0; l<T; l++){
    for(int d=0; d<D; d++){
      distance(d) = abs(theta(l) - theta_dag(d));
    }
    double wl_obs = exponential_work_likeli(theta(l), theta_MLE, J_ni, k);
    int closest = distance.index_min();
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      if(wl_sim(closest,sim) <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec lognormal_pl_mid_P(arma::mat const& theta, arma::mat const& theta_dag, arma::mat const& theta_MLE, int const& k, 
                                 arma::field<arma::mat> const& J_ni, arma::vec const& n_i, int const& M){
  double T = theta.n_rows;
  double D = theta_dag.n_rows;
  arma::mat distance = zeros<mat>(T,D);
  arma::mat wl_sim = zeros<mat>(D,M);
  arma::vec pl = zeros<vec>(T);
  
  for(int d=0; d<D; d++){
    for(int t=0; t<T; t++){
      distance(t,d) = pow(arma::as_scalar(sum(pow(theta.row(t) - theta_dag.row(d),2.0) )), 0.5);
    }
    
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,2);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = R::rlnorm(theta_dag(d,0), pow(theta_dag(d,1), 0.5)); 
        }
        
        theta_MLE_sim(i,0) = mean(log(Ysim_i));
        theta_MLE_sim(i,1) = mean(pow(log(Ysim_i) - theta_MLE_sim(i,0), 2.0));
        J_ni_sim(i) = lognormal_FI(Ysim_i, n_i(i), theta_MLE_sim.row(i).t()); 
      }
      wl_sim(d,sim) = lognormal_work_likeli(theta_dag.row(d).t(), theta_MLE_sim, J_ni_sim, k);
    }
  }
  
  for(int l=0; l<T; l++){
    double wl_obs = lognormal_work_likeli(theta.row(l).t(), theta_MLE, J_ni, k);
    int closest = distance.row(l).index_min();
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      if(wl_sim(closest,sim) <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec logistic_pl_mid_P(arma::mat const& theta, arma::mat const& theta_dag, arma::mat const& theta_MLE, int const& k, 
                             arma::field<arma::mat> const& J_ni, arma::vec const& n_i, 
                             int const& M){
  int p = theta.n_cols;
  double T = theta.n_rows;
  double D = theta_dag.n_rows;
  arma::mat distance = zeros<mat>(T,D);
  arma::mat wl_sim = zeros<mat>(D,M);
  arma::vec pl = zeros<vec>(T);
  
  for(int d=0; d<D; d++){
    for(int t=0; t<T; t++){
      distance(t,d) = pow(arma::as_scalar(sum(pow(theta.row(t) - theta_dag.row(d),2.0) )), 0.5);
    }
    
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,2);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::mat Q, R;
        arma::mat covariates = zeros<mat>(n_i(i),p);
        for(int j=0; j<n_i(i); j++){
          covariates(j,0) = 1;
          covariates(j,1) = R::rnorm(0,1);
        }
        arma::mat Xtheta_dag = covariates * theta_dag.t();
        arma::mat mu_dag = exp(Xtheta_dag)/(1+exp(Xtheta_dag)); 
        arma::qr_econ(Q, R, covariates);
        
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = R::rbinom(1, mu_dag(j,d));  
        }
        
        // GLM optimization modified from package here: github.com/dirkschumacher/rcppglm
        arma::mat C, C2;
        arma::colvec s = arma::zeros<arma::colvec>(p);
        arma::colvec s1, s_old, mu, mu_p, z, W;
        arma::colvec eta = covariates * theta_dag.row(d).t();
        bool is_converged, error;
        error = false;
        
        for (int it = 0; it < 10000; it++) {
          s_old = s;
          mu = (arma::exp(eta) / (1.0 + arma::exp(eta)));
          mu_p = arma::exp(eta) / arma::square(arma::exp(eta) + 1.0);
          z = eta + (Ysim_i - mu) / mu_p;
          W = arma::square(mu_p) / (mu % (1.0 - mu));
          C2 = Q.t() * (Q.each_col() % W);
          if(C2.is_symmetric(1e-16)){
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
          double mu_ij = arma::as_scalar(exp(covariates.row(j)*theta_MLE_sim_i));
          J_ni_sim(i) = J_ni_sim(i) + covariates.row(j).t()*covariates.row(j) * mu_ij / pow(1+mu_ij, 2.0);
        } 
      }
      wl_sim(d,sim) = logistic_work_likeli(theta_dag.row(d).t(), theta_MLE_sim, J_ni_sim, p, k);
    }
  }
  
  for(int l=0; l<T; l++){
    double wl_obs = logistic_work_likeli(theta.row(l).t(), theta_MLE, J_ni, p, k);
    int closest = distance.row(l).index_min();
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      if(wl_sim(closest,sim) <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}

// [[Rcpp::export]]
arma::vec poisson_pl_mid_profile_P(arma::vec const& theta, arma::mat const& theta_dag, double const& meta_estimate, int const& k, 
                                   double const& s_J, arma::vec const& n_i, int const& q, int const& M){
  int p = theta_dag.n_cols;
  double T = theta.size();
  double D = theta_dag.n_rows;
  arma::vec distance = zeros<vec>(D);
  arma::mat wl_sim = zeros<mat>(D,M);
  arma::vec pl = zeros<vec>(T);
  
  for(int d=0; d<D; d++){
    for(int sim=0; sim<M; sim++){
      arma::mat theta_MLE_sim = zeros<mat>(k,p);
      arma::field<arma::mat> J_ni_sim(k);
      for(int i=0; i<k; i++){
        arma::mat Q, R;
        arma::mat covariates = zeros<mat>(n_i(i),p);
        for(int j=0; j<n_i(i); j++){
          covariates(j,0) = 1;
          covariates(j,1) = R::rnorm(0,1);
          covariates(j,2) = R::rnorm(0,1);
        }
        arma::mat Xtheta_dag = covariates * theta_dag.t();
        arma::mat mu_dag = exp(Xtheta_dag); 
        arma::qr_econ(Q, R, covariates);
        
        arma::vec Ysim_i = zeros<vec>(n_i(i));
        
        for(int j=0; j<n_i(i); j++){
          Ysim_i(j) = R::rpois(mu_dag(j,d));  
        }
        
        // GLM optimization modified from package here: github.com/dirkschumacher/rcppglm
        arma::mat C, C2;
        arma::colvec s = arma::zeros<arma::colvec>(p);
        arma::colvec s1, s_old, mu, z;
        arma::colvec eta = covariates * theta_dag.row(d).t();
        bool is_converged;
        
        for (int it = 0; it < 10000; it++) {
          s_old = s;
          mu = arma::exp(eta);
          z = eta + (Ysim_i - mu) / mu;
          C2 = Q.t() * (Q.each_col() % mu);
          if(!C2.is_symmetric(1e-16)) {
            eta = zeros<vec>(n_i(i));
            break;
          }
          try {
            C = arma::chol(C2); 
          } catch(std::runtime_error & e){
            eta = zeros<vec>(n_i(i));
            break;
          }
          s1 = arma::solve(arma::trimatl(C.t()), Q.t() * (mu % z));
          s = arma::solve(arma::trimatu(C), s1);
          eta = Q * s;
          is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < 0.000001;
          if (is_converged) break;
        }
        
        arma::vec theta_MLE_sim_i = arma::solve(arma::trimatu(R), Q.t() * eta);
        theta_MLE_sim.row(i) = theta_MLE_sim_i.t();
        J_ni_sim(i) = zeros<mat>(p,p);
        for(int j=0; j<n_i(i); j++){
          double mu_ij = arma::as_scalar(exp(covariates.row(j)*theta_MLE_sim_i));
          J_ni_sim(i) = J_ni_sim(i) + covariates.row(j).t()*covariates.row(j) * mu_ij;
        } 
      }
      arma::mat s_J_sim = zeros<mat>(p,p);
      arma::vec main_sim = zeros<vec>(p);
      for(int i=0; i<k; i++){
        s_J_sim = s_J_sim + J_ni_sim(i);
        main_sim += J_ni_sim(i) * theta_MLE_sim.row(i).t();
      }
      arma::vec meta_estimate_sim = inv(s_J_sim) * main_sim;
      
      wl_sim(d,sim) = poisson_profile_work_likeli(theta_dag(d,q-1), meta_estimate_sim(q-1), s_J_sim(q-1,q-1));
    }
  }
  
  for(int l=0; l<T; l++){
    for(int d=0; d<D; d++){
      distance(d) = abs(theta(l) - theta_dag(d,q-1));
    }
    double wl_obs = poisson_profile_work_likeli(theta(l), meta_estimate, s_J);
    int closest = distance.index_min();
    arma::vec res = zeros<vec>(M);
    for(int sim=0; sim<M; sim++){
      if(wl_sim(closest,sim) <= wl_obs) res(sim) = 1;
    }
    pl(l) = mean(res);
  }
  return(pl);
}
