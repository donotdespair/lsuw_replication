#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/math/special_functions/bessel.hpp> 
// [[Rcpp::depends(BH)]] 

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
long double omega_conditional_prior (
  double  omega,
  double  prior_s_  = 0.2,
  double  prior_a_  = 1.0
) {
  if (prior_a_ <= 0) stop("'prior_a_' must be greater than 0");
  if (prior_s_ <= 0) stop("'prior_s_' must be greater than 0");
  
  double numerator    = 0.0;
  double denominator  = 0.0;
  double out          = 0.0;
  
  if (omega!= 0) {
    double x      = sqrt(2 / prior_s_) * abs(omega);
    double lambda = prior_a_ - 0.5;
    
    numerator     = pow(abs(omega), prior_a_ - 0.5); // bessel here
    numerator    *= boost::math::cyl_bessel_k(lambda, x);
    denominator   = sqrt(datum::pi);
    denominator  *= pow(2, 0.5 * (prior_a_ - 1.5));
    denominator  *= R::gammafn(prior_a_);
    denominator  *= pow(prior_s_, 0.5 * (prior_a_ + 0.5));
    
    out           = numerator/denominator;
  } else {
    if (prior_a_ <= 0.5){
      out         = datum::inf;
    } else {
      numerator   = R::gammafn(prior_a_ + 1.5);
      denominator = R::gammafn(prior_a_);
      denominator *= sqrt(2 * datum::pi * prior_s_);
      denominator *= pow(prior_a_, 2) - 0.25;
      
      out         = numerator/denominator;
    }
  }
  return out;
} // END omega_conditional_prior



// [[Rcpp::export]]
vec omega_marginal_prior(
    vec     grid_omega, 
    int     S         = 10000,
    double  prior_a_  = 1.0,
    double  prior_s_  = 0.2,
    bool    sample_s_ = true
) {
  
  int G     = grid_omega.n_elem;
  vec marginal_prior(G);
  
  if ( prior_a_ <= 0.5 ) {
    Rcout << "Message: for 'prior$sv_a_' less or equal to 0.5 the ordinate at 0 is unbounded from above." << endl;
  }
  
  if ( sample_s_ ) {
    mat marginal_prior_tmp(G, S);
    vec sample_prior_s_     = prior_s_/as<vec>(Rcpp::rchisq(S, 3));
    
    for (int g=0; g<G; g++) {
      for (int s=0; s<S; s++) {
        marginal_prior_tmp.col(s).row(g) = omega_conditional_prior ( grid_omega(g), sample_prior_s_(s), prior_a_ );
      }
      marginal_prior(g)     = mean(marginal_prior_tmp.row(g));
    }
    
  } else {
    for (int g=0; g<G; g++) {
      marginal_prior(g)     = omega_conditional_prior ( grid_omega(g), prior_s_, prior_a_ );
    }
  }
  
  return marginal_prior;
} // END omega_marginal_prior



// [[Rcpp::export]]
mat  omega_marginal_posterior (
    vec         grid_omega, 
    List&       posterior,  // a list of posteriors
    mat&        Y,          // NxT dependent variables
    mat&        X,          // KxT explanatory variables
    bool        sample_s_ = true
) {
  
  // read inputs
  const cube    posterior_B     = as<cube>(posterior["B"]);
  const cube    posterior_A     = as<cube>(posterior["A"]);
  const cube    posterior_h     = as<cube>(posterior["h"]);
  const cube    posterior_S     = as<cube>(posterior["S"]);
  const mat     posterior_sigma2_omega  = as<mat>(posterior["sigma2_omega"]);
  const mat     posterior_s_    = as<mat>(posterior["s_"]);
  
  const int     S               = posterior_sigma2_omega.n_cols;
  const int     T               = Y.n_cols;
  const int     N               = Y.n_rows;
  const int     G               = grid_omega.n_elem;
  
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  
  mat           marginal_posterior(G, N);
  cube          marginal_posterior_tmp(G, N, S);
  
  for (int g = 0; g < G; g++) {
    for (int n = 0; n < N; n++) {
      for (int s = 0; s < S; s++) {
        mat     residuals       = log(square(posterior_B.slice(s) * (Y - posterior_A.slice(s) * X)));
        
        rowvec  alpha_S(T);
        vec     sigma_S_inv(T);
        for (int t = 0; t < T; t++) {
          rowvec  post_S        = posterior_S.slice(s).row(n);
          alpha_S.col(t)        = alpha_s(post_S(t));
          sigma_S_inv.row(t)    = 1/sigma_s(post_S(t));
        } // END t loop
        
        double  V_omega         = as_scalar(pow(as_scalar(posterior_h.slice(s).row(n) * diagmat(sigma_S_inv) * trans(posterior_h.slice(s).row(n))) + pow(posterior_sigma2_omega(n, s), -1), -1));
        double  omega_bar       = V_omega * as_scalar(posterior_h.slice(s).row(n) * diagmat(sigma_S_inv) * trans(residuals.row(n) - alpha_S));
        marginal_posterior_tmp(g, n, s)   = R::dnorm(grid_omega(g), omega_bar, sqrt(V_omega), false);
      } // END s loop
      
      marginal_posterior(g, n)  = accu(marginal_posterior_tmp.tube(g, n))/S;
    } // END n loop
  } // END s loop 
    
  return(marginal_posterior);
} // END omega_marginal_posterior


/*** R
# omega_conditional_prior(-11, .1, .6)
# ss  = seq(from = -20, to = 20, by = .01)
# ocp = rep(NA, length(ss))
# for (i in 1:length(ss)) ocp[i] = omega_conditional_prior(ss[i], 1, 10)
# plot.ts(ocp, ylim = c(0, max(ocp[-which(ss==0)])))
# omp = omega_marginal_prior(seq(from = -10, to = 10, by = .1))
# plot.ts(omp)

# load("spartan/results/tax_LSUW_2021.rda")
# source("spartan/source_code/thin_bsvar_sv.R")
# qq = thin_sf(qq_BP, 10)
# range(qq_BP$omega)
# ss  = seq(from = -1.1, to = 1.1, by = .01)
# ompo = omega_marginal_posterior(ss, qq, prior, Y, X)
# plot.ts(ompo)
# beepr::beep(4)
*/
