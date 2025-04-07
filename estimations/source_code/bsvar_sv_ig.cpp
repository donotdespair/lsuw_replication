#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppTN)]]
#include "RcppTN.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;

  
/*______________________function orthogonal_complement_matrix_TW______________________*/
mat orthogonal_complement_matrix_TW (const mat& x) {
  // # x is a mxn matrix and m>n
  // # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  // # t(x)%*%out = 0 and det(cbind(x,out))!=0
  int n_nrow     = x.n_rows;
  int n_ncol     = x.n_cols;
  mat Q;
  mat R;
  qr(Q, R, x);
  mat ocm = Q.tail_cols(n_nrow-n_ncol);
  return ocm;
} // END orthogonal_complement_matrix_TW


  
  
/*______________________function do_rgig1______________________*/
// utility function copied from package factorstochvol
double do_rgig1(
    double lambda, 
    double chi, 
    double psi
) { 
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  return as<double>(fun(1, lambda, chi, psi));
} // END do_rgig1



/*______________________function cholesky_tridiagonal______________________*/
// utility function from file precision_sampler.cpp
List cholesky_tridiagonal(
    const vec&    omega_diag,
    const double& omega_offdiag
) {
  const int T = omega_diag.n_elem - 1;
  vec chol_diag(T+1);
  vec chol_offdiag(T+1);
  chol_diag[0] = sqrt(omega_diag[0]);
  for (int j = 1; j < T+1; j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
  return List::create(_["chol_diag"]=chol_diag, _["chol_offdiag"]=chol_offdiag);
} // END cholesky_tridiagonal



/*______________________function forward_algorithm______________________*/
// utility function from file precision_sampler.cpp
vec forward_algorithm(
    const vec& chol_diag,
    const vec& chol_offdiag,
    const vec& covector
) {
  const int T = chol_diag.n_elem - 1;
  vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
} // END forward_algorithm



/*______________________function backward_algorithm______________________*/
// utility function from file precision_sampler.cpp
vec backward_algorithm(
    const vec& chol_diag,
    const vec& chol_offdiag,
    const vec& htmp
) {
  const int T = chol_diag.size() - 1;
  vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
} // END backward_algorithm



/*______________________function precision_sampler_ar1______________________*/
// utility function from file precision_sampler.cpp
vec precision_sampler_ar1(
    const vec&     precision_diag,
    const double&  precision_offdiag,
    const vec&     location
) {
  int T               = location.n_rows;
  vec  epsilon(T, fill::randn);                                                     // sample normal draws 
  List precision_chol = cholesky_tridiagonal(precision_diag, precision_offdiag);    // Cholesky decomposition using a dedicated technique
  vec  aa             = forward_algorithm(precision_chol["chol_diag"],              // this forward substitution can be performed outside of the loop
                                          precision_chol["chol_offdiag"],
                                                        location);
  vec draw_ssar1      = backward_algorithm(precision_chol["chol_diag"],
                                           precision_chol["chol_offdiag"],
                                                         aa + epsilon);     // this has to be done in the loop as function backward_algorithm requires covector to be a vector (not a matrix)
  return draw_ssar1;
} // END precision_sampler_ar1



/*______________________function inverse_transform_sampling______________________*/
// utility function from file utils_latent_states.cc from the source code of package stochvol
uvec inverse_transform_sampling (
    const vec&  mixprob,
    const int   T
) {
  uvec r(T);
  for (int j = 0; j < T; j++) {
    int index = (10-1)/2;  // start searching in the middle
    const double unnorm_cdf_value = R::unif_rand()*mixprob[9 + 10*j];  // current (non-normalized) value
    bool larger = false;  // indicates that we already went up
    bool smaller = false; // indicates that we already went down
    while(true) {
      if (unnorm_cdf_value > mixprob[index +  10*j]) {
        index++;
        if (smaller) {
          break;
        } else {
          larger = true;
        }
      } else if (larger || index == 0) {
        break;
      } else {
        index--;
        smaller = true;
      }
    }
    r[j] = index;
  }
  return r;
}



/*______________________function find_mixture_indicator_cdf______________________*/
// utility function from file utils_latent_states.cc from the source code of package stochvol
vec find_mixture_indicator_cdf (
    const vec& datanorm           // provide all that is conditionally normal
){
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  
  const int T = datanorm.n_elem;
  vec mixprob(10 * T);
  for (int j = 0; j < T; j++) {  // TODO slow (10*T calls to exp)!
    const int first_index = 10*j;
    mixprob(first_index) = std::exp(pr_s(0) - (datanorm(j) - alpha_s(0)) * (datanorm(j) - alpha_s(0)) / sigma_s(0) );
    for (int r = 1; r < 10; r++) {
      mixprob(first_index+r) = mixprob(first_index+r-1) + std::exp(pr_s(r) - (datanorm(j) - alpha_s(r)) * (datanorm(j) - alpha_s(r)) / sigma_s(r) );
    }
  }
  return mixprob;
}



/*______________________function svar_ce1______________________*/
List svar_ce1 (
    rowvec&         aux_h_n,            // 1xT
    double&         aux_rho_n,
    double&         aux_omega_n,
    double&         aux_sigma2v_n,
    double&         aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&         aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    urowvec&        aux_S_n,            // 1xT
    const rowvec&   u,                  // 1xT
    const List&     prior,
    bool            sample_s_ = true
) {
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const rowvec  U = log(pow(u + ccc, 2));
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(trans(U - aux_omega_n*aux_h_n));
  aux_S_n               = trans(inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // // sample aux_s_n
  // if ( sample_s_ ) {
  //   aux_s_n               = (1 + 2 * aux_sigma2_omega_n) / chi2rnd(3 + 2 * prior_sv_a_);
  // }
  
  // // sample aux_sigma2_omega
  // aux_sigma2_omega_n    = randg( distr_param(1 + 0.5 * prior_sv_a_, pow(pow(prior_sv_s_,-1) + pow(2 * aux_sigma2v_n,-1), -1)  ) );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow( hm1 * hm1.t() / aux_sigma2v_n, -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * (hm1 * aux_h_n.cols(1,T-1).t() / aux_sigma2v_n) );
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-1,1);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_sigma2v
  aux_sigma2v_n         = (aux_sigma2_omega_n + as_scalar(aux_h_n * HH_rho * aux_h_n.t())) / chi2rnd( 3 + T );
  aux_omega_n           = pow(aux_sigma2v_n, 0.5);
  
  // sample aux_h
  mat       V_h         = diagmat(sigma_S_inv) + (HH_rho / aux_sigma2v_n);
  vec       h_bar       = diagmat(sigma_S_inv) * (U - alpha_S).t();
  aux_h_n               = trans(precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END svar_ce1



/*______________________function sample_A_heterosk1 ______________________*/
void sample_A_heterosk1 (
    mat&        aux_A,          // NxK
    const mat&  aux_B,          // NxN
    const vec&  aux_hyper,      // NxM
    const mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const mat&  Y,              // NxT dependent variables
    const mat&  X,              // KxT dependent variables
    const List& prior           // a list of priors - original dimensions
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = pow(aux_hyper(1), -1) * as<mat>(prior["A_V_inv"]);
  rowvec    zerosA(K);
  vec sigma_vectorised= vectorise(aux_sigma);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   zn_sigma    = zn / sigma_vectorised;
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
    
    mat     precision = prior_A_Vinv + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * prior_A_Vinv + trans(zn_sigma) * Wn_sigma;
    
    mat     precision_chol = trimatu(chol(precision));
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + as<vec>(rnorm(K)));
    aux_A.row(n)      = trans(draw);
  } // END n loop
} // END sample_A_heterosk1




/*______________________function sample_B_heterosk1______________________*/
void sample_B_heterosk1 (
    mat&        aux_B,          // NxN
    const mat&  aux_A,          // NxK
    const vec&  aux_hyper,      // NxM
    const mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const mat&  Y,              // NxT dependent variables
    const mat&  X,              // KxT dependent variables
    const List& prior,          // a list of priors - original dimensions
    const field<mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  mat prior_SS_inv          = pow(aux_hyper(0), -1) * as<mat>(prior["B_V_inv"]);
  mat shocks                = Y - aux_A * X;
  
  
  for (int n=0; n<N; n++) {
    
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    mat posterior_SS_inv    = prior_SS_inv + shocks_sigma * shocks_sigma.t();
    mat posterior_S_inv     = VB(n) * posterior_SS_inv * VB(n).t();
    posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
    
    // sample B
    mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u                 = rnorm(posterior_nu+1, 0, pow(posterior_nu, -0.5));
    alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn                = Rcpp::rnorm(rn-1, 0, pow(posterior_nu, -0.5));
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)           = b0n * VB(n);
  } // END n loop
} // END sample_B_heterosk1



/*______________________function sample_hyperparameters______________________*/
void sample_hyperparameters (
    vec&              aux_hyper,
    const mat&        aux_B,
    const mat&        aux_A,
    const field<mat>& VB,
    const List&       prior,
    LogicalVector     sample_hyper
) {
  // the function changes the value of aux_hyper by reference (filling it with a new draw)
  const int N = aux_B.n_rows;
  const int K = aux_A.n_cols;
  int rn=0;
  for (int n=0; n<N; n++) {
    rn       += VB(n).n_rows;
  }
  
  if (sample_hyper[4]) {
    aux_hyper(4)    = ( as<double>(prior["hyper_S"]) + aux_hyper(2) + aux_hyper(3) ) / R::rchisq(as<double>(prior["hyper_V"]) + 4*as<double>(prior["hyper_a"]));
  }
  if (sample_hyper[3]) {
    aux_hyper(3)    = R::rgamma( as<double>(prior["hyper_a"]) + 0.5*as<double>(prior["hyper_nu"]) ,
              1/((1/aux_hyper(4)) + (1/(2*aux_hyper(1)))));
  }
  if (sample_hyper[2]) {
    aux_hyper(2)    = R::rgamma( as<double>(prior["hyper_a"]) + 0.5*as<double>(prior["hyper_nu"]) ,
              1/((1/aux_hyper(4)) + (1/(2*aux_hyper(0)))) );
  }
  if (sample_hyper[1]) {
    aux_hyper(1)    = ( aux_hyper(3) + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
      R::rchisq( as<double>(prior["hyper_nu"]) + N * K );
  }
  if (sample_hyper[0]) {
    aux_hyper(0)    = ( aux_hyper(2) + trace(aux_B * as<mat>(prior["B_V_inv"]) * trans(aux_B) )) /
      R::rchisq( as<double>(prior["hyper_nu"]) + rn );
  }
} // END sample_hyperparameters



  
  
  
// [[Rcpp::export]]
Rcpp::List bsvar_svce_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,         // restrictions on B0
    const Rcpp::List&             starting_values, 
    const bool                    sample_s_ = true,
    const int                     thin = 10,  // introduce thinning
    const LogicalVector           sample_hyper = LogicalVector::create(true, true, false, false, false)
) {
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR-SV model              |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
  Rcout << "    Every " << thin << "th draw is saved via MCMC thinning" << endl;
  Rcout << " Press Esc to interrupt the computations" << endl;
  Rcout << "**************************************************|" << endl;
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  
  const int   S     = floor(SS / thin);
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  vec   aux_hyper   = as<vec>(starting_values["hyper"]);  // 5x1 (gamma_0, gamma_+, s_0, s_+, s_)
  mat   aux_h       = as<mat>(starting_values["h"]);
  vec   aux_rho     = as<vec>(starting_values["rho"]);
  vec   aux_omega   = as<vec>(starting_values["omega"]);
  vec   aux_sigma2v = as<vec>(starting_values["sigma2v"]);
  umat  aux_S       = as<umat>(starting_values["S"]);
  vec   aux_sigma2_omega = as<vec>(starting_values["sigma2_omega"]);
  vec   aux_s_      = as<vec>(starting_values["s_"]);
  mat   aux_sigma(N, T);
  
  for (int n=0; n<N; n++) {
    aux_sigma.row(n) = exp(0.5 * aux_omega(n) * aux_h.row(n));
  }
  
  cube  posterior_B(N, N, S);
  cube  posterior_A(N, K, S);
  mat   posterior_hyper(5, S);
  cube  posterior_h(N, T, S);
  mat   posterior_rho(N, S);
  mat   posterior_omega(N, S);
  mat   posterior_sigma2v(N, S);
  ucube posterior_S(N, T, S);
  mat   posterior_sigma2_omega(N, S);
  mat   posterior_s_(N, S);
  cube  posterior_sigma(N, T, S);
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_hyper
    sample_hyperparameters( aux_hyper, aux_B, aux_A, VB, prior, sample_hyper);
    
    // sample aux_B
    sample_B_heterosk1(aux_B, aux_A, aux_hyper, aux_sigma, Y, X, prior, VB);
    
    // sample aux_A
    sample_A_heterosk1(aux_A, aux_B, aux_hyper, aux_sigma, Y, X, prior);
    
    // sample aux_h, aux_omega and aux_S, aux_sigma2_omega
    mat U = aux_B * (Y - aux_A * X);
    
    for (int n=0; n<N; n++) {
      rowvec  h_tmp     = aux_h.row(n);
      double  rho_tmp   = aux_rho(n);
      double  omega_tmp = aux_omega(n);
      double  sigma2v_tmp = aux_sigma2v(n);
      urowvec S_tmp     = aux_S.row(n);
      rowvec  U_tmp     = U.row(n);
      double  s2o_tmp   = aux_sigma2_omega(n);
      double  s_n       = aux_s_(n);
      
      List sv_n         = svar_ce1( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, U_tmp, prior, sample_s_ );
      
      aux_h.row(n)      = as<rowvec>(sv_n["aux_h_n"]);
      aux_rho(n)        = as<double>(sv_n["aux_rho_n"]);
      aux_omega(n)      = as<double>(sv_n["aux_omega_n"]);
      aux_sigma2v(n)    = as<double>(sv_n["aux_sigma2v_n"]);
      aux_S.row(n)      = as<urowvec>(sv_n["aux_S_n"]);
      aux_sigma2_omega(n) = as<double>(sv_n["aux_sigma2_omega_n"]);
      aux_s_(n)         = as<double>(sv_n["aux_s_n"]);
      
      aux_sigma.row(n)  = exp(0.5 * aux_h.row(n));
    }
    
    if (ss % thin == 0) {
      posterior_B.slice(s)          = aux_B;
      posterior_A.slice(s)          = aux_A;
      posterior_hyper.col(s)        = aux_hyper;
      posterior_h.slice(s)          = aux_h;
      posterior_rho.col(s)          = aux_rho;
      posterior_omega.col(s)        = aux_omega;
      posterior_sigma2v.col(s)      = aux_sigma2v;
      posterior_S.slice(s)          = aux_S;
      posterior_sigma2_omega.col(s) = aux_sigma2_omega;
      posterior_s_.col(s)           = aux_s_;
      posterior_sigma.slice(s)      = aux_sigma;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["h"]        = aux_h,
      _["rho"]      = aux_rho,
      _["omega"]    = aux_omega,
      _["sigma2v"]  = aux_sigma2v,
      _["S"]        = aux_S,
      _["sigma2_omega"] = aux_sigma2_omega,
      _["s_"]       = aux_s_,
      _["sigma"]    = aux_sigma
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["h"]        = posterior_h,
      _["rho"]      = posterior_rho,
      _["omega"]    = posterior_omega,
      _["sigma2v"]  = posterior_sigma2v,
      _["S"]        = posterior_S,
      _["sigma2_omega"] = posterior_sigma2_omega,
      _["s_"]        = posterior_s_,
      _["sigma"]    = posterior_sigma
    )
  );
} // END bsvar_svce_cpp
