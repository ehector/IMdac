// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector Cpp_Qgk_log_deriv(NumericVector z, double A, double B, double g, double k, double c){
  // Some of this is borrowed from the R package "gk"
  int n = z.length();
  NumericVector A_vec(n,A);
  NumericVector B_vec(n,B);
  NumericVector g_vec(n,g);
  NumericVector k_vec(n,k);
  NumericVector c_vec(n,c);
  NumericVector z_squared = pow(z, 2.0);
  NumericVector term1 = k_vec * log(1 + z_squared);
  NumericVector term2 = 1 + c_vec * tanh(g_vec * z/2);
  NumericVector term3 = (1 + (2 * k_vec + 1) * z_squared)/(1 + z_squared);
  NumericVector term4 = c_vec * g_vec * z/(2 * pow(cosh(g_vec * z/2), 2.0));
  LogicalVector gzero = (g_vec == 0);
  for(int q=0; q<n; q++){
    if(gzero(q)) {
      term2(q) = 1;
      term4(q) = 0;
    }
    if(std::isnan(z_squared(q))){
      term1(q) = 2 * k * log(std::abs(z(q)));
      term3(q) = 2 * k + 1;
      term4(q) = 0;
    }
    if(k_vec(q)==0) term1(q) = 0;
  }
  return(log(B_vec) + term1 + log(term2 * term3 + term4));
}

// [[Rcpp::export]]
NumericVector g_and_k_MLE(NumericVector Ysim_i,  double a0,  double A,  double alpha,  NumericVector c0,  double gamma,  int p, 
                          NumericVector theta0,  NumericVector theta_min,  NumericVector theta_max, int  maxit,  double tol){
  // Some of this is borrowed from the R package "gk"
  
  NumericVector theta_est = theta0;
  NumericVector theta_prev(p, 0.0);
  
  Rcpp::Function pgk = Environment::global_env()["pgk"];
  
  double diff_ll, hatL1, hatL2;
  bool convergence = false;
  const bool log_flag = TRUE;
  int it = 0;
  while(it<maxit+1 && !convergence){
    double at = a0 * pow(it + 1 + A, -alpha);
    NumericVector ct = c0 * pow(it + 1, -gamma);
    NumericVector gt(p, 0.0);
    
    for(int r=0; r<p; r++){
      NumericVector delta(p,0.0);
      delta(r) = 1.0;
      NumericVector theta1 = pmin(pmax(theta_est + ct*delta, theta_min), theta_max);
      NumericVector theta2 = pmin(pmax(theta_est - ct*delta, theta_min), theta_max);
      
      SEXP z_SEXP1 = pgk(Ysim_i, theta1(0), theta1(1), theta1(2), theta1(3), 0.8, log_flag);
      NumericVector z_vec1 = as<NumericVector>(z_SEXP1);
      NumericVector z_norm1 = Rcpp::dnorm(z_vec1, 0, 1, log_flag);
      SEXP z_SEXP2 = pgk(Ysim_i, theta2(0), theta2(1), theta2(2), theta2(3), 0.8, log_flag);
      NumericVector z_vec2 = as<NumericVector>(z_SEXP2);
      NumericVector z_norm2 = Rcpp::dnorm(z_vec2, 0, 1, log_flag);
      NumericVector hatL1_vec = z_norm1 - Cpp_Qgk_log_deriv(z_vec1, theta1(0), theta1(1), theta1(2), theta1(3), 0.8);
      NumericVector hatL2_vec = z_norm2 - Cpp_Qgk_log_deriv(z_vec2, theta2(0), theta2(1), theta2(2), theta2(3), 0.8);
      
      hatL1 = -sum(hatL1_vec);
      hatL2 = -sum(hatL2_vec);
      
      gt(r) = (hatL1 - hatL2)/(theta1(r) - theta2(r));
    }
    theta_est = pmin(pmax(theta_est - at * gt, theta_min), theta_max);
    
    if(it==0) {
      for(int r=0; r<p; r++){
        theta_prev(r) = theta_est(r);
      }
      it += 1;
    } else {
      diff_ll= sum(abs(theta_est - theta_prev));
      if(diff_ll > tol) {
        for(int r=0; r<p; r++){
          theta_prev(r) = theta_est(r);
        }
        it += 1;
      } else {
        convergence=true;
      }
    }
  }
  return(theta_est);
}

