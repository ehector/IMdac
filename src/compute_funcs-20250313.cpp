// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double exponential_work_likeli(double const& theta, arma::vec const& theta_MLE, arma::vec const& J_ni, int const& k){
  double s_J = sum(J_ni);
  double theta_n = arma::as_scalar(J_ni.t() * theta_MLE)/s_J;
  double ll = exp(arma::as_scalar( - s_J * pow(theta_n - theta, 2.0)) / 2);
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
