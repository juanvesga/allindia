
#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::vec compute_dx(const arma::vec& invec, 
                     arma::mat& lambda, 
                     arma::uvec lo_id,
                     arma::uvec hi_id,
                     const arma::mat& lin, 
                     const arma::mat& lo_ds, 
                     const arma::mat& lo_mdr, 
                     const arma::mat& hi_ds,
                     const arma::mat& hi_mdr,
                     const arma::vec& mortvec,
                     const arma::mat& smear, 
                     const arma::mat& xpert,
                     const arma::mat& xray,
                     const arma::mat& acf,
                     double pgrowth, 
                     arma::uvec birth_id,
                     arma::uvec nores,
                     arma::uvec respSymp,
                     double pr_rs, 
                     double pr_hi, 
                     const arma::mat& agg_inc,
                     const arma::mat& agg_notif,
                     const arma::mat& agg_mort,
                     const arma::mat& agg_pt,
                     const arma::mat& agg_txfp,
                     const arma::mat& agg_dx,
                     const arma::mat& agg_acf,
                     const arma::mat& agg_pmo,
                     const arma::mat& sel_inc,
                     const arma::mat& sel_acqu,
                     const arma::mat& sel_remo,
                     const arma::mat& sel_notif,
                     const arma::mat& sel_pt,
                     const arma::mat& sel_txfp,
                     const arma::mat& sel_dx,
                     const arma::mat& sel_acf,
                     arma::uvec aux_inc,
                     arma::uvec aux_notif,
                     arma::uvec aux_mort,
                     arma::uvec aux_remo,
                     arma::uvec aux_pt,
                     arma::uvec aux_txfp,
                     arma::uvec aux_dx,
                     arma::uvec aux_pmo,
                     arma::uvec aux_acf,
                     int n_aux,
                     double t) {
  
  arma::vec lam = lambda*invec;
  
  // -- Get all model components together
  arma::mat allmat (lin +  
    lam(0) * lo_ds  +
    lam(1) * lo_mdr +
    lam(2) * hi_ds  +
    lam(3) * hi_mdr );
  
  
  arma::vec dx = allmat*invec;
  
  // Implement deaths
  arma::vec morts = mortvec%invec;
  dx = dx - morts;
  //  Implement births
  double births = sum(morts(nores))*(pgrowth==0) + pgrowth;
  double births_rs = sum(morts(respSymp))*(pgrowth==0) + pgrowth*pr_rs;
  
  dx(birth_id(0)) = dx(birth_id(0))+births*(1-pr_hi) ;
  dx(birth_id(1)) = dx(birth_id(1))+births*pr_hi ;
  dx(birth_id(2)) = dx(birth_id(2))+births_rs*(1-pr_hi) ;
  dx(birth_id(3)) = dx(birth_id(3))+births_rs*pr_hi ;
  
  
  // Get the auxiliaries
  arma::vec ou(n_aux, arma::fill::zeros);
  arma::vec out;  
  out= arma::join_cols(dx,ou);
  out(aux_inc)         = agg_inc*(sel_inc%allmat)*invec;
  out(aux_inc(2))      = out(aux_inc(2)) + arma::accu((sel_acqu%allmat)*invec);
  out(aux_remo)        = agg_inc.rows(0,2)*(sel_remo%lin)*invec;
  out(aux_notif)       = agg_notif*(sel_notif%allmat)*invec;
  out(aux_mort)        = agg_mort*morts;
  
  // Get auxuliares for counting units (only for intervention times)
  if (t>=2017){
 
    out(aux_pt)         = agg_pt*(sel_pt%allmat)*invec;
    out(aux_txfp)       = agg_txfp*(sel_txfp%allmat)*invec;
    out(aux_dx(0))      = arma::accu(agg_dx*(sel_dx%smear)*invec);
    out(aux_dx(1))      = arma::accu(agg_dx*(sel_dx%xpert)*invec);
    out(aux_dx(2))      = arma::accu(agg_dx*(sel_dx%xray)*invec);
    out(aux_pmo)        = agg_pmo*invec;
    out(aux_acf)        = agg_acf*(sel_acf%acf)*invec;
    
    
  }
  
  return out;
  
}

