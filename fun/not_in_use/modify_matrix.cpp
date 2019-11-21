
#include <RcppArmadillo.h>


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::sp_mat modify_matrix(const arma::sp_mat& lin, const arma::sp_mat& linh, 
                           const arma::vec& lam,
                  const arma::sp_mat& a0ds, const arma::sp_mat& a0mdr,
                  const arma::sp_mat& a5ds, const arma::sp_mat& a5mdr,
                  const arma::sp_mat& a15ds, const arma::sp_mat& a15mdr,
                  const arma::sp_mat& a65ds, const arma::sp_mat& a65mdr,
                 const arma::sp_mat& nlinh, double foi) {


// -- Get all model components together
  arma::sp_mat M (lin + linh +
    lam(0) * a0ds  +
    lam(1) * a0mdr +
    lam(2) * a5ds  +
    lam(3) * a5mdr +
    lam(4) * a15ds  +
    lam(5) * a15mdr +
    lam(6) * a65ds  +
    lam(7) * a65mdr +
    foi * nlinh);
  
  return M;

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be autosp_matically 
// run after the compilation.
//

