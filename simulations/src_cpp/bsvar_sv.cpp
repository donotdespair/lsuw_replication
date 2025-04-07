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
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat orthogonal_complement_matrix_TW (const arma::mat& x) {
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


// [[Rcpp::export]]
std::string ordinal(
    int n
) {
  std::string suffix;
  if (n % 10 == 1 && n % 100 != 11) {
    suffix = "st";
  } else if (n % 10 == 2 && n % 100 != 12) {
    suffix = "nd";
  } else if (n % 10 == 3 && n % 100 != 13) {
    suffix = "rd";
  } else {
    suffix = "th";
  }
  return std::to_string(n) + suffix;
} // END ordinal

















// [[Rcpp::export]]
double do_rgig1(
    double lambda, 
    double chi, 
    double psi
) { 
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  return as<double>(fun(1, lambda, chi, psi));
} // END do_rgig1


// [[Rcpp::export]]
Rcpp::List cholesky_tridiagonal(
    const arma::vec&    omega_diag,
    const double&       omega_offdiag
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


// [[Rcpp::export]]
arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector
) {
  const int T = chol_diag.n_elem - 1;
  vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
} // END forward_algorithm


// [[Rcpp::export]]
arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp
) {
  const int T = chol_diag.size() - 1;
  vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
} // END backward_algorithm


// [[Rcpp::export]]
arma::vec precision_sampler_ar1(
    const arma::vec&     precision_diag,
    const double&        precision_offdiag,
    const arma::vec&     location
) {
  int T               = location.n_rows;
  vec  epsilon(T, fill::randn);
  List precision_chol = cholesky_tridiagonal(precision_diag, precision_offdiag);    // Cholesky decomposition using a dedicated technique
  vec  aa             = forward_algorithm(precision_chol["chol_diag"],              // this forward substitution can be performed outside of the loop
                                          precision_chol["chol_offdiag"],
                                                        location);
  vec draw_ssar1      = backward_algorithm(precision_chol["chol_diag"],
                                           precision_chol["chol_offdiag"],
                                                         aa + epsilon);     // this has to be done in the loop as function backward_algorithm requires covector to be a vector (not a matrix)
  return draw_ssar1;
} // END precision_sampler_ar1


// [[Rcpp::export]]
arma::uvec inverse_transform_sampling (
    const arma::vec&  mixprob,
    const int         T
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



// [[Rcpp::export]]
arma::vec find_mixture_indicator_cdf (
    const arma::vec& datanorm           // provide all that is conditionally normal
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


// [[Rcpp::export]]
Rcpp::List svar_nc1 (
    arma::rowvec&   aux_h_n,            // 1xT
    double&         aux_rho_n,
    double&         aux_omega_n,
    double&         aux_sigma2v_n,
    double&         aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&         aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&  aux_S_n,            // 1xT
    const arma::rowvec&   u,                  // 1xT
    const Rcpp::List&     prior,
    bool            sample_s_ = true
) {
  // sampler for the non-centred parameterisation of the SV process
  
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
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (prior_sv_s_ + 2 * aux_sigma2_omega_n)/chi2rnd(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = do_rgig1( prior_sv_a_-0.5, pow(aux_omega_n,2), 2/aux_s_n );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow(hm1*hm1.t(), -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  double    upper_bound = pow(1-aux_sigma2_omega_n, 0.5);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_omega
  double    V_omega_inv = 1/( as_scalar(aux_h_n * diagmat(sigma_S_inv) * aux_h_n.t()) + pow(aux_sigma2_omega_n, -1) );
  double    omega_bar   = as_scalar(aux_h_n * diagmat(sigma_S_inv) * (U - alpha_S).t());
  double    omega_aux   = randn( distr_param(V_omega_inv*omega_bar, sqrt(V_omega_inv) ));
  
  // sample aux_h
  mat       V_h         = pow(omega_aux, 2) * diagmat(sigma_S_inv) + HH_rho;
  vec       h_bar       = omega_aux * diagmat(sigma_S_inv) * (U - alpha_S).t();
  rowvec    h_aux       = trans(precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  // ASIS
  rowvec    aux_h_tilde = omega_aux * h_aux;
  double    hHHh        = as_scalar( aux_h_tilde * HH_rho * aux_h_tilde.t() );
  aux_sigma2v_n         = do_rgig1( -0.5*(T-1), hHHh, 1/aux_sigma2_omega_n );
  int       ss=1;
  if (R::runif(0,1)<0.5) ss *= -1;
  aux_omega_n           = ss * sqrt(aux_sigma2v_n);
  aux_h_n               = aux_h_tilde / aux_omega_n;
  
  // ASIS: resample aux_rho
  hm1                   = aux_h_n.cols(0,T-2);
  aux_rho_var           = as_scalar(pow(hm1*hm1.t(), -1));
  aux_rho_mean          = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  upper_bound           = pow(1-aux_sigma2_omega_n, 0.5);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END sv_nc1



// [[Rcpp::export]]
Rcpp::List svar_ce1 (
    arma::rowvec&       aux_h_n,            // 1xT
    double&             aux_rho_n,
    double&             aux_omega_n,
    double&             aux_sigma2v_n,
    double&             aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&             aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&      aux_S_n,            // 1xT
    const arma::rowvec& u,                  // 1xT
    const Rcpp::List&   prior,
    bool                sample_s_ = true
) {
  // sampler for the centred parameterisation of the SV process
  
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
  
  // sample aux_s_n
  // if ( sample_s_ ) {
  //   aux_s_n               = (1 + 2 * aux_sigma2_omega_n) / chi2rnd(3 + 2 * prior_sv_a_);
  // }
  
  // sample aux_sigma2_omega
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






















// [[Rcpp::export]]
arma::mat sample_B_heterosk1 (
    arma::mat&        aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat shocks                = Y;
  
  
  for (int n=0; n<N; n++) {
    
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    mat posterior_SS_inv    = (pow(aux_hyper(n,0), -1) * prior_SS_inv) + shocks_sigma * shocks_sigma.t();
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
    vec   u(posterior_nu+1, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)           = b0n * VB(n);
  } // END n loop
  
  return aux_B;
} // END sample_B_heterosk1






// [[Rcpp::export]]
arma::mat sample_hyperparameters (
    arma::mat&              aux_hyper,       // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&        aux_B,            // NxN
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior
) {
  // the function returns aux_hyper by reference (filling it with a new draw)
  
  const int N = aux_B.n_rows;
  // const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  // double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  // double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  // double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  // double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  // mat   prior_A               = as<mat>(prior["A"]);
  // mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  aux_hyper(2 * N, 0) = scale_tmp / R::rchisq(shape_tmp);
  
  // aux_A - related hyper-parameters 
  // ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  // scale_tmp           = prior_hyper_s_AA + 2 * sum(ss_tmp);
  // shape_tmp           = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  // aux_hyper(2 * N, 1) = scale_tmp / R::rchisq(shape_tmp);
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn            = VB(n).n_rows;
    
    // aux_B - related hyper-parameters 
    scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
    shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
    aux_hyper(N + n, 0) = R::rgamma(shape_tmp, scale_tmp);
    
    scale_tmp         = aux_hyper(N + n, 0) + as_scalar(aux_B.row(n) * prior_B_V_inv * aux_B.row(n).t());
    shape_tmp         = prior_hyper_nu_B + rn;
    aux_hyper(n, 0)   = scale_tmp / R::rchisq(shape_tmp);
    
    // // aux_A - related hyper-parameters 
    // scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
    // shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
    // aux_hyper(N + n, 1) = R::rgamma(shape_tmp, scale_tmp);
    // 
    // scale_tmp         = aux_hyper(N + n, 1) + 
    //   as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
    // shape_tmp         = prior_hyper_nu_A + K;
    // aux_hyper(n, 1)   = scale_tmp / R::rchisq(shape_tmp);
  } // END n loop
  
  return aux_hyper;
} // END sample_hyperparameters





// [[Rcpp::export]]
Rcpp::List normaliseB_s (
    mat& aux_B,
    mat& aux_hyper,
    mat& aux_h,
    vec& aux_rho,
    vec& aux_omega,
    vec& aux_sigma2v,
    umat& aux_S,
    vec& aux_sigma2_omega,
    vec& aux_s_,
    mat& aux_sigma
) {
  // only works for 2x2 B matrix
  mat Bsq = square(aux_B);
  
  if (Bsq(0,0) < Bsq(0,1)) {
    // create permutation matrix
    umat P(2, 2);
    P(0,1) = 1;
    P(1,0) = 1;
    
    int N   = aux_B.n_rows;
    // permutation
    aux_B     = P * aux_B;
    aux_hyper.rows(0, N - 1) = P * aux_hyper.rows(0, N - 1);
    aux_hyper.rows(N, 2 * N - 1) = P * aux_hyper.rows(N, 2 * N - 1);
    aux_h     = P * aux_h;
    aux_rho   = P * aux_rho;
    aux_omega = P * aux_omega;
    aux_sigma2v = P * aux_sigma2v;
    aux_S     = P * aux_S;
    aux_sigma2_omega = P * aux_sigma2_omega;
    aux_s_    = P * aux_s_;
    aux_sigma = P * aux_sigma;
  }
  
  return List::create(
    _["aux_B"]  = aux_B,
    _["aux_hyper"] = aux_hyper,
    _["aux_h"]  = aux_h,
    _["aux_rho"] = aux_rho,
    _["aux_omega"] = aux_omega,
    _["aux_sigma2v"] = aux_sigma2v,
    _["aux_S"]  = aux_S,
    _["aux_sigma2_omega"] = aux_sigma2_omega,
    _["aux_s_"] = aux_s_,
    _["aux_sigma"] = aux_sigma
  );
} // END normaliseB_s





// [[Rcpp::export]]
Rcpp::List bsvar_sv_cpp (
    const int&                    S,          // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,         // restrictions on B0
    const Rcpp::List&             starting_values, 
    const int                     thin = 100, // introduce thinning
    const bool                    centred_sv = false,
    const bool                    show_progress = true
) {
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = ordinal(thin) + " ";
  }
  
  std::string       name_model = "";
  if ( centred_sv ) {
    name_model        = "    Centred";
  } else {
    name_model        = "Non-centred";
  }
  
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Gibbs sampler for the SVAR-SV model              |" << endl;
    Rcout << "   " << name_model << " SV model is estimated              |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  // const int   K     = X.n_rows;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  // mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);  
  mat   aux_h       = as<mat>(starting_values["h"]);
  vec   aux_rho     = as<vec>(starting_values["rho"]);
  vec   aux_omega   = as<vec>(starting_values["omega"]);
  vec   aux_sigma2v = as<vec>(starting_values["sigma2v"]);
  umat  aux_S       = as<umat>(starting_values["S"]);
  vec   aux_sigma2_omega = as<vec>(starting_values["sigma2_omega"]);
  vec   aux_s_      = as<vec>(starting_values["s_"]);
  mat   aux_sigma(N, T);
  
  if ( centred_sv ) {
    for (int n=0; n<N; n++) {
      aux_sigma.row(n) = exp(0.5 * aux_h.row(n));
    }
  } else {
    for (int n=0; n<N; n++) {
      aux_sigma.row(n) = exp(0.5 * aux_omega(n) * aux_h.row(n));
    }
  }
  
  const int   SS     = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  // cube  posterior_A(N, K, SS);
  cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  posterior_h(N, T, SS);
  mat   posterior_rho(N, SS);
  mat   posterior_omega(N, SS);
  mat   posterior_sigma2v(N, SS);
  ucube posterior_S(N, T, SS);
  mat   posterior_sigma2_omega(N, SS);
  mat   posterior_s_(N, SS);
  cube  posterior_sigma(N, T, SS);
  
  int   ss = 0;
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_hyper
    aux_hyper       = sample_hyperparameters( aux_hyper, aux_B, VB, prior);
    
    // sample aux_B
    aux_B           = sample_B_heterosk1(aux_B, aux_hyper, aux_sigma, Y, prior, VB);
    
    // sample aux_A
    // aux_A           = sample_A_heterosk1(aux_A, aux_B, aux_hyper, aux_sigma, Y, X, prior);
    
    // sample aux_h, aux_omega and aux_S, aux_sigma2_omega
    // mat U = aux_B * (Y - aux_A * X);
    mat U = aux_B * Y;
    
    for (int n=0; n<N; n++) {
      rowvec  h_tmp     = aux_h.row(n);
      double  rho_tmp   = aux_rho(n);
      double  omega_tmp = aux_omega(n);
      double  sigma2v_tmp = aux_sigma2v(n);
      urowvec S_tmp     = aux_S.row(n);
      rowvec  U_tmp     = U.row(n);
      double  s2o_tmp   = aux_sigma2_omega(n);
      double  s_n       = aux_s_(n);
      
      List sv_n;
      if ( centred_sv ) {
        sv_n            = svar_ce1( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, U_tmp, prior, true );
      } else {
        sv_n            = svar_nc1( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, U_tmp, prior, true );
      }

      aux_h.row(n)      = as<rowvec>(sv_n["aux_h_n"]);
      aux_rho(n)        = as<double>(sv_n["aux_rho_n"]);
      aux_omega(n)      = as<double>(sv_n["aux_omega_n"]);
      aux_sigma2v(n)    = as<double>(sv_n["aux_sigma2v_n"]);
      aux_S.row(n)      = as<urowvec>(sv_n["aux_S_n"]);
      aux_sigma2_omega(n)         = as<double>(sv_n["aux_sigma2_omega_n"]);
      aux_s_(n)         = as<double>(sv_n["aux_s_n"]);

      if ( centred_sv ) {
        aux_sigma.row(n)  = exp(0.5 * aux_h.row(n));
      } else {
        aux_sigma.row(n)  = exp(0.5 * aux_omega(n) * aux_h.row(n));
      }
    }
    
    // List ll             = normaliseB_s ( aux_B, aux_hyper, aux_h, aux_rho, aux_omega, aux_sigma2v, aux_S, aux_sigma2_omega, aux_s_, aux_sigma);
    // 
    // aux_B             = as<mat>(ll["aux_B"]); 
    // aux_hyper         = as<mat>(ll["aux_hyper"]);
    // aux_h             = as<mat>(ll["aux_h"]);
    // aux_rho           = as<vec>(ll["aux_rho"]);
    // aux_omega         = as<vec>(ll["aux_omega"]);
    // aux_sigma2v       = as<vec>(ll["aux_sigma2v"]);
    // aux_S             = as<umat>(ll["aux_S"]);
    // aux_sigma2_omega  = as<vec>(ll["aux_sigma2_omega"]);
    // aux_s_            = as<vec>(ll["aux_s_"]);
    // aux_sigma         = as<mat>(ll["aux_sigma"]);
    
    if (s % thin == 0) {
      posterior_B.slice(ss)          = aux_B;
      // posterior_A.slice(ss)          = aux_A;
      posterior_hyper.slice(ss)      = aux_hyper;
      posterior_h.slice(ss)          = aux_h;
      posterior_rho.col(ss)          = aux_rho;
      posterior_omega.col(ss)        = aux_omega;
      posterior_sigma2v.col(ss)       = aux_sigma2v;
      posterior_S.slice(ss)          = aux_S;
      posterior_sigma2_omega.col(ss) = aux_sigma2_omega;
      posterior_s_.col(ss)           = aux_s_;
      posterior_sigma.slice(ss)      = aux_sigma;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      // _["A"]        = aux_A,
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
      // _["A"]        = posterior_A,
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
} // END bsvar_sv_cpp








// [[Rcpp::export]]
arma::vec log_mean (
    arma::mat     log_density     // n x s matrix with log density ordinates
) {
  int S               = log_density.n_cols;
  vec c_log_density   = max(log_density, 1);
  vec log_numerator   = c_log_density - log(S) + log( sum( exp(log_density.each_col() - c_log_density), 1) );
  return log_numerator;
} // log_mean 





// [[Rcpp::export]]
arma::vec verify_volatility_sv_cpp (
    const Rcpp::List&       posterior,  // a list of posteriors
    const Rcpp::List&       prior,      // a list of priors - original dimensions
    const arma::mat&        Y,          // NxT dependent variables
    const bool              sample_s_ = true
) {
  // computes the log of SDDR for homoskedasticity hypothesis omega_n = 0
  // see Lütkepohl, Shang, Uzeda, Woźniak (2013)
  
  // read inputs
  const cube    posterior_B     = as<cube>(posterior["B"]);
  const cube    posterior_h     = as<cube>(posterior["h"]);
  const cube    posterior_S     = as<cube>(posterior["S"]);
  const mat     posterior_sigma2_omega  = as<mat>(posterior["sigma2_omega"]);
  const mat     posterior_s_    = as<mat>(posterior["s_"]);
  
  const double  prior_a_        = as<double>(prior["sv_a_"]);
  const double  prior_s_        = as<double>(prior["sv_s_"]);
  
  const int     S               = posterior_sigma2_omega.n_cols;
  const int     T               = Y.n_cols;
  // const int     N               = Y.n_rows;
  
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  
  if ( prior_a_ <= 0.5 ) {
    stop("'prior$sv_a_' must be greater than 0.5");
  }
  
  // compute denominator
  double inv_sqrt_s_      = 0.0;
  vec sample_prior_s_(S);
  if ( sample_s_ ) {
    sample_prior_s_       = prior_s_/chi2rnd( 3, S );
    inv_sqrt_s_           = as_scalar(mean(pow(sample_prior_s_, -0.5)));
  } else {
    inv_sqrt_s_           = pow(prior_s_, -0.5);
  }
  double  log_denominator     = - 0.5 * log(2 * M_PI) + log(inv_sqrt_s_) - log(pow(prior_a_, 2) - 0.25) + R::lgammafn(prior_a_ + 1.5) - R::lgammafn(prior_a_);
  
  // compute numerator
  mat     log_numerator_s(1, S);
  for (int s = 0; s < S; s++) {
    int n = 0;
    // for (int n = 0; n < N; n++) {
    mat     residuals       = log(square(posterior_B.slice(s) * Y ));
      
      rowvec  alpha_S(T);
      vec     sigma_S_inv(T);
      
      for (int t = 0; t < T; t++) {
        rowvec  post_S        = posterior_S.slice(s).row(n);
        alpha_S.col(t)        = alpha_s(post_S(t));
        sigma_S_inv.row(t)    = 1/sigma_s(post_S(t));
      } // END t loop
      
      double  V_omega         = pow(as_scalar(posterior_h.slice(s).row(n) * diagmat(sigma_S_inv) * trans(posterior_h.slice(s).row(n))) + pow(posterior_sigma2_omega(n, s), -1), -1);
      double  omega_bar       = V_omega * as_scalar(posterior_h.slice(s).row(n) * diagmat(sigma_S_inv) * trans(residuals.row(n) - alpha_S));
      log_numerator_s(n, s)   = R::dnorm(0, omega_bar, sqrt(V_omega), true);
    // } // END n loop
  } // END s loop
  
  // compute the log of the mean numerator exp(log_numerator)
  vec log_numerator           = log_mean(log_numerator_s);

  return log_numerator - log_denominator;
} // END verify_volatility_cpp


// [[Rcpp::export]]
double normalise_jaro_distance (
    mat&          B,              // NxN
    const mat&    B_hat,          // NxN
    const mat&    sigma2,         // NxT
    const int     T,
    const mat&    communication   // NNxNN communication matrix
) {
  const int N     = B.n_rows; 
  vec   vecBt     = vectorise(B.t());
  vec   vecB_hatt = vectorise(B_hat.t());
  
  mat   B_hat_inv = inv(B_hat);
  vec   sigma2_sum = sum(sigma2, 1);
  
  mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
  mat   covariance_inv  = T * Finsher_info;
  
  double distance = as_scalar(trans(vecBt - vecB_hatt) * covariance_inv * (vecBt - vecB_hatt));
  
  return distance;
} // END normalise_jaro_distance



// [[Rcpp::export]]
mat covariance_inv_bsvars (
    const mat&    B_hat,          // NxN
    const mat&    sigma2,         // NxT
    const int     T,
    const mat&    communication   // NNxNN communication matrix
) {
  const int N     = B_hat.n_rows; 
  vec   vecB_hatt = vectorise(B_hat.t());
  mat   B_hat_inv = inv(B_hat);
  vec   sigma2_sum = sum(sigma2, 1);
  
  mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
  mat   covariance_inv  = T * Finsher_info;
  
  return covariance_inv;
} // END covariance_inv_bsvars


// [[Rcpp::export]]
mat pick_permutation_bsvars (
    arma::mat&              B,
    const arma::mat&        B_hat,
    const int               T,              // sample size
    const arma::mat&        cov_inv,
    const arma::umat&       perms,
    const arma::mat&        signs,
    const arma::mat&        communication   // NNxNN communication matrix
) {
  
  const int   no_perms    = perms.n_rows;
  const int   no_signs    = signs.n_rows;
  const int   N           = B.n_rows;
  
  mat   I                 = eye(N,N);
  vec   vecB_hatt         = vectorise(B_hat.t());
  vec   vecBt(vecB_hatt.n_elem);
  
  mat   distance(no_perms, no_signs);
  
  for (int p=0; p<no_perms; p++) {
    for (int s=0; s<no_signs; s++) {
      vecBt               = vectorise(trans(I.rows(perms.row(p)) * diagmat(signs.row(s)) * B));
      distance(p, s)      = as_scalar(trans(vecBt - vecB_hatt) * cov_inv * (vecBt - vecB_hatt));
    } // END s loop
  } // END p loop
  uvec sign_min_inds      = index_min(distance, 1);
  vec sign_mins           = min(distance, 1);
  int perm_min_ind        = index_min(sign_mins);
  int sign_min_ind        = sign_min_inds(perm_min_ind);
  
  mat P = I.rows(perms.row(perm_min_ind));
  mat D = diagmat(trans(signs.row(sign_min_ind)));
  
  return P * D;
  // return distance;
} // END pick_permutation_bsvars






// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_sv (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  // cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  cube   posterior_h       = as<cube>(posterior["h"]);
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  mat    posterior_rho     = as<mat>(posterior["rho"]);
  mat    posterior_omega   = as<mat>(posterior["omega"]);
  mat    posterior_sigma2v = as<mat>(posterior["sigma2v"]);
  cube   posterior_S       = as<cube>(posterior["S"]);
  mat    posterior_sigma2_omega = as<mat>(posterior["sigma2_omega"]);
  mat    posterior_s_      = as<mat>(posterior["s_"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  Function AllPerms("allPermutations");
  Function AllSigns("allSigns");
  Function communication_mat("communication_matrix");
  
  mat communication = as<mat>(communication_mat(N, N));
  
  mat cov_inv(pow(N,2), pow(N,2));
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    mat sigma2      = square(posterior_sigma.slice(s));
    cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_h.slice(s)      = Pabs * posterior_h.slice(s);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_rho.col(s)      = Pabs * posterior_rho.col(s);
    posterior_omega.col(s)    = Pabs * posterior_omega.col(s);
    posterior_sigma2v.col(s)  = Pabs * posterior_sigma2v.col(s);
    posterior_S.slice(s)      = Pabs * posterior_S.slice(s);
    posterior_sigma2_omega.col(s) = Pabs * posterior_sigma2_omega.col(s);
    posterior_s_.col(s)       = Pabs * posterior_s_.col(s);
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
  } // END s loop
  
  return List::create(
    _["B"]            = posterior_B,
    _["hyper"]        = posterior_hyper,
    _["h"]            = posterior_h,
    _["sigma"]        = posterior_sigma,
    _["rho"]          = posterior_rho,
    _["omega"]        = posterior_omega,
    _["sigma2v"]      = posterior_sigma2v,
    _["S"]            = posterior_S,
    _["sigma2_omega"] = posterior_sigma2_omega,
    _["s_"]           = posterior_s_
  );
} // END normalise_jaro_bsvar_sv


/*** R


allPermutations = function(N) {
  return(rbind(1:N, permute::allPerms(N)))
}

allSigns = function(N) {
  diag.signs    = matrix(NA, 2^N, N)
  for (n in 1:N) {
    diag.signs[,n]    = kronecker(c(-1,1), rep(1, 2^(n - 1)))
  }
  return(diag.signs)
}

communication_matrix = function(N, M) {
  stopifnot("N must be a positive integer." = N %% 1 == 0 & N > 0)
  stopifnot("M must be a positive integer." = M %% 1 == 0 & M > 0)
  
  NM = N * M
  
  return(diag(NM)[as.vector(t(matrix(1:NM, N, M))),])
}


# options(echo=TRUE)
# args <- commandArgs(trailingOnly = TRUE)
# args
# iteration <- as.integer(args[1])
# iteration
# rm(args)
# 
# iteration = 1
# library(bsvars)
# set.seed(123 + iteration)
# S_burn          = 1e3
# S               = 1e3
# N               = 2
# T               = 260
# yy              = matrix(rnorm(N * T), T, N)
# 
# B               = matrix(TRUE, N, N)
# spec            = specify_bsvar_sv$new(yy, B = B)
# prior           = spec$prior$get_prior()
# starting_values = spec$starting_values$get_starting_values()
# starting_values$h = matrix(rnorm(prod(dim(yy)), sd = 0.01), ncol(yy), nrow(yy))
# starting_values$S = matrix(1, ncol(yy), nrow(yy))
# starting_values$B = matrix(c(10, -10, 0, 20), N, N)
# starting_values$omega = c(0, sqrt(0.05))
# starting_values$rho = c(0, 0.9)
# VB              = spec$identification$get_identification()
# 
# burn            = bsvar_sv_cpp(S = S_burn, Y = t(yy), prior, VB, starting_values, thin = 1, centred_sv = FALSE, show_progress = TRUE)
# post            = bsvar_sv_cpp(S = S, Y = t(yy), prior, VB, starting_values = burn$last_draw, thin = 1, centred_sv = FALSE, show_progress = TRUE)
# par(mfrow = c(2, 2)); for (i in 1:2) { for (j in 1:2) {hist(post$posterior$B[i,j,],breaks=100)}}
# # B_hat = matrix(c(10, -10, 0, 20), N, N)
# B_hat = diag(diag(sign(post$last_draw$B))) %*% post$last_draw$B
# post$posterior = normalise_jaro_bsvar_sv (post$posterior, B_hat = B_hat)
# par(mfrow = c(2, 2)); for (i in 1:2) { for (j in 1:2) {hist(post$posterior$B[i,j,],breaks=100)}}
# apply(post$posterior$B, 1:2, mean)
# sddr            = -as.numeric(verify_volatility_sv_cpp( post$posterior, prior, t(yy) ))
# save(sddr, file = paste0("bsvar_sv_",iteration,".rda"))

*/
