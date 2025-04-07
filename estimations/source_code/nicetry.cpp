#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat nicetry (
    irowvec  sign,
    urowvec  perm,
    mat   A
) {
  
  urowvec  perm_cpp  = perm - 1;
  mat   I         = eye(3,3);
  
  return I.rows(perm_cpp) * diagmat(trans(sign)) * A;
}
    

/*** R
sig = c(-1,1,1)
perm = c(2,1,3)
A = matrix(1:9, 3, 3)
nicetry(sig,perm, A)

*/
