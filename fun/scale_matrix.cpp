
#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::mat scale_matrix(const arma::mat& m0, const arma::mat& m1, 
                       double scalefac) {
  
  
  // -- Get all model components together
  arma::mat M (m0 + scalefac*(m1-m0) );
  
  return M;
  
}


