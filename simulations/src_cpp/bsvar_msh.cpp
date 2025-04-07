
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;



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
arma::vec log_mean (
    arma::mat     log_density     // n x s matrix with log density ordinates
) {
  int S               = log_density.n_cols;
  vec c_log_density   = max(log_density, 1);
  vec log_numerator   = c_log_density - log(S) + log( sum( exp(log_density.each_col() - c_log_density), 1) );
  return log_numerator;
} // log_mean 




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
int csample_num1 (
    Rcpp::NumericVector x,
    Rcpp::NumericVector prob = NumericVector::create()
) {
  bool replace = false;
  NumericVector ret = Rcpp::RcppArmadillo::sample(x, 1, replace, prob);
  int out           = ret(0);
  return out;
} // END csample_num1





// [[Rcpp::export]]
arma::mat count_regime_transitions (
    const arma::mat& xi
) {
  const int M = xi.n_rows;
  const int T = xi.n_cols;
  
  mat count(M, M);
  urowvec s   = index_max( xi, 0 );
  for (int t=1; t<T; t++) {
    count( s(t-1), s(t))++;
  }
  return count;
} // END count_regime_transitions



// [[Rcpp::export]]
arma::rowvec rDirichlet1 (
    const arma::rowvec&   alpha     // Kx1
) {
  const int K   = alpha.size();
  rowvec    draw(K);
  for (int k=0; k<K; k++) {
    draw(k)     = randg(distr_param(alpha(k), 1.0));
  }
  return draw/sum(draw);
} // END rDirichlet1



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
Rcpp::List sample_transition_probabilities (
    arma::mat           aux_PR_TR,    // MxM 
    arma::vec           aux_pi_0,     // Mx1
    const arma::mat&    aux_xi,       // MxT
    const Rcpp::List&   prior,         // a list of priors - original dimensions
    const bool          MSnotMIX = true
) {
  // the function changes the value of aux_PR_TR and aux_pi_0 by reference (filling it with a new draw)
  const int   M           = aux_PR_TR.n_rows;
  const mat   prior_PR_TR = as<mat>(prior["PR_TR"]);
  
  if ( MSnotMIX ) {
    mat transitions       = count_regime_transitions(aux_xi);
    mat posterior_alpha   = transitions + prior_PR_TR;
    
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = rDirichlet1(posterior_alpha.row(m));
    }
    vec prob_xi1          = aux_PR_TR *aux_xi.col(0);
    prob_xi1             /= sum(prob_xi1);
    int S0_draw           = csample_num1(wrap(seq_len(M)), wrap(prob_xi1));
    rowvec posterior_alpha_0(M, fill::value((1.0)));
    posterior_alpha_0(S0_draw-1)++;
    aux_pi_0              = trans(rDirichlet1(posterior_alpha_0));
  } else {
    rowvec occurrences    = trans(sum(aux_xi, 1));
    rowvec posterior_alpha= occurrences + prior_PR_TR.row(0);
    aux_pi_0              = trans(rDirichlet1(posterior_alpha));
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = aux_pi_0.t();
    }
  }
  
  return List::create(
    _["PR_TR"]        = aux_PR_TR,
    _["pi_0"]         = aux_pi_0
  );
} // END sample_transition_probabilities









// [[Rcpp::export]]
arma::mat filtering_msh (
    const arma::mat&  U,                  // NxT
    const arma::mat&  sigma,              // NxM
    const arma::mat&  PR_TR,              // MxM
    const arma::vec&  pi_0                // Mx1
) {
  const int   T = U.n_cols;
  const int   N = U.n_rows;
  const int   M = PR_TR.n_rows;
  
  mat         eta_t(M, T);
  mat         xi_t_t(M, T);
  
  // This loop evaluates mvnormal pdf at U - simplified operations for zero-mean diagonal-covariance case
  for (int m=0; m<M; m++) {
    rowvec log_d    = -0.5 * sum(pow( pow(sigma.col(m), -0.5) % U.each_col(), 2), 0);
    log_d          += -0.5 * N * log(2*M_PI) - 0.5 * log(prod(sigma.col(m)));
    NumericVector   exp_log_d   = wrap(exp(log_d));
    exp_log_d[exp_log_d==0]     = 1e-300;
    eta_t.row(m)    = as<rowvec>(exp_log_d);
  } // END m loop
  
  vec xi_tm1_tm1    = pi_0;
  
  for (int t=0; t<T; t++) {
    vec     num     = eta_t.col(t) % (PR_TR.t() * xi_tm1_tm1);
    double  den     = sum(num);
    xi_t_t.col(t)   = num/den;
    xi_tm1_tm1      = xi_t_t.col(t);
  } // END t loop
  
  return xi_t_t;
} // END filtering_msh




// [[Rcpp::export]]
arma::mat smoothing_msh (
    const arma::mat&  U,                  // NxT
    const arma::mat&  PR_TR,              // MxM
    const arma::mat&  filtered            // MxT
) {
  const int   T = U.n_cols;
  const int   M = PR_TR.n_rows;
  
  mat   smoothed(M, T);
  smoothed.col(T-1)   = filtered.col(T-1);
  
  for (int t=T-2; t>=0; --t) {
    smoothed.col(t)   = (PR_TR * (smoothed.col(t+1)/(PR_TR.t() * filtered.col(t)) )) % filtered.col(t);
    if (any(smoothed.col(t) < 0) || any(smoothed.col(t) > 1)) {
      for (int m=0; m<M; m++) {
        if (smoothed(m,t) > 1) {smoothed(m,t) = 1;}
        if (smoothed(m,t) < 0) {smoothed(m,t) = 0;}
      }
      smoothed.col(t) = smoothed.col(t) / accu(smoothed.col(t));
    }
  } // END t loop
  
  return smoothed;
} // smoothing_msh




// [[Rcpp::export]]
arma::mat sample_Markov_process_msh (
    arma::mat&        aux_xi,             // MxT
    const arma::mat&  U,                  // NxT
    const arma::mat&  aux_sigma2,         // NxM
    const arma::mat&  aux_PR_TR,         // MxM
    const arma::vec&  aux_pi_0,          // Mx1
    const bool        finiteM = true
) {
  // the function changes the value of aux_xi by reference (filling it with a new draw)
  
  int minimum_regime_occurrences = 0;
  int max_iterations = 1;
  if ( finiteM ) {
    minimum_regime_occurrences = 2;
    max_iterations = 10;
  }
  
  const int   T   = U.n_cols;
  const int   M   = aux_PR_TR.n_rows;
  mat aux_xi_tmp  = aux_xi;
  
  mat filtered    = filtering_msh(U, aux_sigma2, aux_PR_TR, aux_pi_0);
  mat smoothed    = smoothing_msh(U, aux_PR_TR, filtered);
  mat     aj      = eye(M, M);
  
  mat xi(M, T);
  int draw        = csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
  aux_xi_tmp.col(T-1)     = aj.col(draw-1);
  
  if ( minimum_regime_occurrences==0 ) {
    for (int t=T-2; t>=0; --t) {
      vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
      draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
      aux_xi_tmp.col(t)   = aj.col(draw-1);
    }
    aux_xi = aux_xi_tmp;
  } else {
    int regime_occurrences  = 1;
    int iterations  = 1;
    while ( (regime_occurrences<minimum_regime_occurrences) & (iterations<max_iterations) ) {
      for (int t=T-2; t>=0; --t) {
        vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
        draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
        aux_xi_tmp.col(t)   = aj.col(draw-1);
      }
      mat transitions       = count_regime_transitions(aux_xi_tmp);
      regime_occurrences    = min(transitions.diag());
      iterations++;
    } // END while
    if ( iterations<max_iterations ) aux_xi = aux_xi_tmp;
  }
  
  return aux_xi;
} // END sample_Markov_process_msh





// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::mat sample_variances_msh (
    arma::mat&          aux_sigma2, // NxM
    const arma::mat&    aux_B,      // NxN
    const arma::mat&    Y,          // NxT dependent variables
    const arma::mat&    aux_xi,     // MxT state variables
    const Rcpp::List&   prior       // a list of priors - original dimensions
) {
  // the function changes the value of aux_sigma2 by reference (filling it with a new draw)
  const int   M     = aux_xi.n_rows;
  const int   N     = aux_B.n_rows;
  const int   T     = Y.n_cols;
  
  rowvec posterior_nu   = sum(aux_xi, 1).t() + as<double>(prior["sigma_nu"]);
  mat posterior_s(N, M);
  posterior_s.fill(prior["sigma_s"]);
  for (int m=1; m<M; m++) {
    
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        posterior_s.col(m) += square(aux_B * Y.col(t));
      }
    }
    
    for (int n=0; n<N; n++) {
      aux_sigma2(n, m)   = posterior_s(n, m) /  chi2rnd(posterior_nu(n));
    }
    
  }
  
  return aux_sigma2;
} // END sample_variances_msh






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

  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  mat    prior_B_V_inv        = as<mat>(prior["B_V_inv"]);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  aux_hyper(2 * N, 0) = scale_tmp / R::rchisq(shape_tmp);
  
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
    
  } // END n loop
  
  return aux_hyper;
} // END sample_hyperparameters






// [[Rcpp::export]]
void normaliseB_s (
    mat& aux_B,
    mat& aux_hyper,
    mat& aux_sigma2, 
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
    aux_sigma2 = P * aux_sigma2;
    aux_sigma = P * aux_sigma;
  }
} // END normaliseB_s




// [[Rcpp::export]]
Rcpp::List bsvar_msh_cpp (
    const int&              S,              // No. of posterior draws
    const arma::mat&        Y,              // NxT dependent variables
    const Rcpp::List&       prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,       // restrictions on B0
    const Rcpp::List&       starting_values,
    const int               thin = 100,     // introduce thinning
    const bool              finiteM = true,
    const bool              show_progress = true
) {
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = ordinal(thin) + " ";
  }
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
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
  mat   aux_sigma2  = as<mat>(starting_values["sigma2"]);
  mat   aux_sigma(N, T);
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  
  const int   M     = aux_PR_TR.n_rows;
  
  const int   SS     = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  // cube  posterior_A(N, K, SS);
  cube  posterior_sigma2(N, M, SS);
  cube  posterior_PR_TR(M, M, SS);
  mat   posterior_pi_0(M, SS);
  cube  posterior_xi(M, T, SS);
  cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  posterior_sigma(N, T, SS);
  
  int   ss = 0;
  for (int t=0; t<T; t++) {
    aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi.col(t).index_max()) , 0.5 );
  }
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_hyper DONE
    aux_hyper         = sample_hyperparameters(aux_hyper, aux_B, VB, prior);
    
    // sample aux_B DONE
    aux_B             = sample_B_heterosk1(aux_B, aux_hyper, aux_sigma, Y, prior, VB);
    
    // // sample aux_A DONE
    // aux_A             = sample_A_heterosk1(aux_A, aux_B, aux_hyper, aux_sigma, Y, X, prior);
    
    // sample aux_xi
    mat U = aux_B * Y;
    aux_xi            = sample_Markov_process_msh(aux_xi, U, aux_sigma2, aux_PR_TR, aux_pi_0, finiteM);
    
    // sample aux_PR_TR DONE
    List aux_PR_tmp   = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior, true);
    aux_PR_TR         = as<mat>(aux_PR_tmp["PR_TR"]);
    aux_pi_0          = as<vec>(aux_PR_tmp["pi_0"]);
    
    // // sample aux_sigma2
    aux_sigma2        = sample_variances_msh(aux_sigma2, aux_B, Y, aux_xi, prior);
    
    for (int t=0; t<T; t++) {
      aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi.col(t).index_max()) , 0.5 );
    }
    
    if (s % thin == 0) {
      posterior_B.slice(ss)      = aux_B;
      // posterior_A.slice(ss)      = aux_A;
      posterior_sigma2.slice(ss) = aux_sigma2;
      posterior_PR_TR.slice(ss)  = aux_PR_TR;
      posterior_pi_0.col(ss)     = aux_pi_0;
      posterior_xi.slice(ss)     = aux_xi;
      posterior_hyper.slice(ss)  = aux_hyper;
      posterior_sigma.slice(ss)  = aux_sigma;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      // _["A"]        = aux_A,
      _["sigma2"]   = aux_sigma2,
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi,
      _["hyper"]    = aux_hyper,
      _["sigma"]    = aux_sigma
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      // _["A"]        = posterior_A,
      _["sigma2"]   = posterior_sigma2,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi,
      _["hyper"]    = posterior_hyper,
      _["sigma"]    = posterior_sigma
    )
  );
} // END bsvar_msh



// [[Rcpp::export]]
double dig2r (
    const double x, 
    const double a1, 
    const double a2, 
    const double b1, 
    const double b2,
    const bool   logarithm = false
) {
  
  double lc = 0.5 * a1 * log( b1 );
  lc       += 0.5 * a2 * log( b2 );
  lc       -= R::lbeta(0.5 * a1, 0.5 * a2);
  
  double lk = 0.5 * (a2 - 2) * log(x);
  lk       -= 0.5 * (a1 + a2) * log(b1 + b2 * x);
  
  double out = lc + lk;
  if ( !logarithm ) {
    out = exp( out );
  }
  return out;
} // END dig2r




// [[Rcpp::export]]
arma::vec verify_volatility_msh_cpp (
    const Rcpp::List&       posterior,  // a list of posteriors
    const Rcpp::List&       prior,      // a list of priors - original dimensions
    const arma::mat&        Y          // NxT dependent variables
) {
  // code working only for N=2 and M=2
  
  cube  posterior_sigma2    = posterior["sigma2"];
  cube  posterior_B         = posterior["B"];
  cube  posterior_xi        = posterior["xi"];
  
  const int   M             = posterior_xi.n_rows;
  const int   N             = posterior_B.n_rows;
  const int   T             = Y.n_cols;
  const int   S             = posterior_B.n_slices;
  
  rowvec      homoskedasticity_hypothesis(M, fill::ones);
  
  // compute denominator
  
  double  sigma_nu          = as<double>(prior["sigma_nu"]);
  double  sigma_s           = as<double>(prior["sigma_s"]);
  
  double  log_denominator = dig2r( 1, sigma_nu, sigma_nu, sigma_s, sigma_s, true);
  
  // compute numerator
  vec     log_numerator_s(S);
  for (int s = 0; s < S; s++) {
      
      rowvec posterior_nu   = sum(posterior_xi.slice(s), 1).t() + sigma_nu;
      
      mat posterior_s(N, M);
      posterior_s.fill(prior["sigma_s"]);
      for (int m=0; m<M; m++) {
        for (int t=0; t<T; t++) {
          if (posterior_xi(m,t,s)==1) {
            posterior_s.col(m) += square(posterior_B.slice(s) * Y.col(t) );
          }
        }
      }
      
      log_numerator_s(s)  = dig2r( 1, posterior_nu(0), posterior_nu(1), posterior_s(0), posterior_s(1), true );
      
  } // END s loop
  
  // compute the log of the mean numerator exp(log_numerator)
  vec log_numerator           = log_mean(trans(log_numerator_s));
  
  return log_numerator - log_denominator;
} // END verify_volatility_msh_cpp




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
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_sigma2  = as<cube>(posterior["sigma2"]);
  cube   posterior_PR_TR   = as<cube>(posterior["PR_TR"]);
  cube   posterior_xi      = as<cube>(posterior["xi"]);
  mat    posterior_pi_0    = as<mat>(posterior["pi_0"]);
  
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
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_sigma2.slice(s) = Pabs * posterior_sigma2.slice(s);
  
  } // END s loop
  
  return List::create(
    _["B"]            = posterior_B,
    _["hyper"]        = posterior_hyper,
    _["sigma"]        = posterior_sigma,
    _["sigma2"]       = posterior_sigma2,
    _["PR_TR"]        = posterior_PR_TR,
    _["pi_0"]         = posterior_pi_0,
    _["xi"]           = posterior_xi
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

*/



/*** R
# 
# options(echo=TRUE)
# args <- commandArgs(trailingOnly = TRUE)
# iteration <- as.integer(args[1])
# rm(args)
# 
# library(bsvars)
# set.seed(123 + iteration)
# S_burn          = 1e1
# S               = 1e2
# N               = 2
# T               = 260
# yy              = matrix(rnorm(N * T), T, N)
# 
# M               = 2
# B               = matrix(TRUE, N, N)
# spec            = specify_bsvar_msh$new(yy, B = B, p = 1, M = M)
# prior           = spec$prior$get_prior()
# prior$PR_TR     = matrix(1, M, M) + 9 * diag(M)
# starting_values = spec$starting_values$get_starting_values()
# starting_values$xi = cbind(starting_values$xi[,1], starting_values$xi)
# VB              = spec$identification$get_identification()
# 
# burn  = bsvar_msh_cpp ( S_burn, t(yy), prior, VB, starting_values, thin = 1)
# burn$posterior$B[1,,]^2
# post  = bsvar_msh_cpp ( S, t(yy), prior, VB, burn$last_draw, thin = 1)
# 
# 
# sddr            = -as.numeric(verify_volatility_msh_cpp( post$posterior, prior, t(yy) ))
# save(sddr, file = paste0("bsvar_msh_",iteration,".rda"))

*/
