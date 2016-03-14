/** momErr.cpp
 * 
 * Code to compute the first and second moments of the errors of a simulation
 * Philip Barrett
 * Chicago, 09mar2016
 */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <string.h>
#include <boost/math/distributions/normal.hpp>
#include <math.h>

using namespace Rcpp ;
using namespace arma ;


// [[Rcpp::export]]
arma::vec ar1_sim( int n_pds, double rho, double sig_eps, 
            bool init_flag=false, double init=0 ){
// Creates an AR(1) simulation
  vec innov = sig_eps * randn( n_pds ) ;
      //  The innovation vector
  vec out = zeros(n_pds) ;
      //  Initialize the output
  if( init_flag ){
    out(0) = init ;
        // Initialize the output
  }
  else
  {
    double sig_uncond = sig_eps / std::sqrt( 1 - std::pow( rho, 2.0 ) ) ;
        // The unconditional variance
    vec v_init = sig_uncond * randn( 1 ) ;
    out(0) = v_init(0) ;
        // Fill the first element of the output with a random draw from the
        // unconditional distribution
  }
  for( int i = 1 ; i < n_pds ; i++ ){
    out(i) = rho * out(i-1) + innov(i) ;
        // Fill the output with the AR(1) process
  }
  return out ;
}

// [[Rcpp::export]]
arma::vec norm_thresh( int n_pds, double rho, double sig_eps ){
// Produces draws from a random normal with the unconditional variance of an
// AR(1) process
  double sig_uncond = sig_eps / std::sqrt( 1 - std::pow( rho, 2.0 ) ) ;
  vec out = sig_uncond * randn( n_pds ) ;
  return out ;
}

// [[Rcpp::export]]
arma::ivec signal_create( arma::vec theta, arma::vec theta_hat ){
// Returns a matrix of signals
  int n_pds = theta.n_elem ;
  ivec out = zeros<ivec>( n_pds ) ;
  for( int i = 0 ; i < n_pds ; i++ ){
    out(i) = ( theta(i) > theta_hat(i) ) ? 1 : 0 ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::rowvec mu_sig2_update( arma::rowvec mu_sig2, double theta_hat,
                              double sig_eps, double rho, int y ){
// Computes the updated mu, sigma^2 pair
  double mu = mu_sig2(0) ;
  double sig2 = mu_sig2(1) ;
  double sig = sqrt( sig2 ) ;
  double psi = ( theta_hat - mu ) / sig ;
      // Basics
  boost::math::normal norm_dist( 0, 1 ) ;
      // Instantiate a standard normal
  double num = pdf( norm_dist, psi ) ;
      // The numerator
  double denom = ( y == 1 ) ? 1 - cdf( norm_dist, psi ) : cdf( norm_dist, psi ) ;
      // The denominator
  double haz_raz ;
      // The (reverse) hazard rate
  double sign = ( y == 1 ) ? 1 : -1 ;
  if( denom > 1e-13 ){
    haz_raz = num / denom ;
  }else{
    // Use the asymptotic expansion if the denominator is too small
        // Upper or lower tail
    haz_raz = sign * ( psi + 1 / psi - 2 / pow( psi, 3 ) + 10 / pow( psi, 5 ) 
                          - 74 / pow( psi, 7 ) + 706 / pow( psi, 9 ) ) ;
  }
  double mu_update = rho * ( mu + sign * sig * haz_raz) ;
  double sig2_update = pow( sig_eps, 2 ) + pow( rho, 2 ) * sig2 *
                        ( 1 - pow( haz_raz, 2 ) + sign * psi * haz_raz ) ;
      // The updated distribution parameters
  rowvec out(2) ;
  out << mu_update << sig2_update ;
      // Format the output
  return out ;
}

// [[Rcpp::export]]
arma::mat thresh_filter( arma::rowvec mu_sig2_0, arma::rowvec theta_hat, 
                            double sig_eps, double rho, arma::irowvec y ){
// Apply the threshold filter
  int K = y.n_elem ;
      // Number of periods
  mat out = zeros( K + 1, 2 ) ;
  out.row(0) = mu_sig2_0 ;
      // Initialize the output matrix
  for( int i=0; i < K ; i++ ){
    out.row(i+1) = mu_sig2_update( out.row(i), theta_hat(i), sig_eps, rho, y(i) ) ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::rowvec mu_sig2_update_gf( arma::rowvec mu_sig2, double theta_hat,
                              double sig_eps, double rho, int y ){
// Computes the updated mu, sigma^2 pair for the Gaussian filter
  double mu = mu_sig2(0) ;
  double sig2 = mu_sig2(1) ;
  double sig = sqrt( sig2 ) ;
  double psi = ( theta_hat - mu ) / sig ;
      // Basics
  boost::math::normal norm_dist( 0, 1 ) ;
      // Instantiate a standard normal
  double num = pdf( norm_dist, psi ) ;
      // The numerator
  double denom_h = 1 - cdf( norm_dist, psi ) ;
  double denom_r = cdf( norm_dist, psi ) ;
      // The denominators
  double haz = num / denom_h ;
  double raz = num / denom_r ;
      // The (reverse) hazard rate
  double sign = ( y == 1 ) ? 1 : -1 ;
  if( denom_h < 1e-13 ){
    haz = psi + 1 / psi - 2 / pow( psi, 3 ) + 10 / pow( psi, 5 ) 
                          - 74 / pow( psi, 7 ) + 706 / pow( psi, 9 ) ;
  }     // Numerical correction for the hazard rate
  if( denom_r < 1e-13 ){
    raz = - ( psi + 1 / psi - 2 / pow( psi, 3 ) + 10 / pow( psi, 5 ) 
                          - 74 / pow( psi, 7 ) + 706 / pow( psi, 9 ) ) ;
  }     // Numerical correction for the reverse hazard rate
  double mu_update = rho * ( mu + sig * ( ( y == 1 ) ? haz : - raz ) ) ;
  double sig2_update = pow( sig_eps, 2 ) + pow( rho, 2 ) * sig2 * 
                              ( 1 - haz * raz ) ;
      // The updated distribution parameters
  rowvec out(2) ;
  out << mu_update << sig2_update ;
      // Format the output
  return out ;
}


// [[Rcpp::export]]
arma::mat gauss_filter( arma::rowvec mu_sig2_0, arma::rowvec theta_hat, 
                            double sig_eps, double rho, arma::irowvec y ){
// Apply the threshold filter
  int K = y.n_elem ;
      // Number of periods
  mat out = zeros( K + 1, 2 ) ;
  out.row(0) = mu_sig2_0 ;
      // Initialize the output matrix
  for( int i=0; i < K ; i++ ){
    out.row(i+1) = mu_sig2_update_gf( out.row(i), theta_hat(i), sig_eps, rho, y(i) ) ;
  }
  return out ;
}


// [[Rcpp::export]]
arma::mat multi_ar1_sim( int n_sim, int n_pds, double rho, double x0, 
                          double sig_eps ){
// Code to create multiple simulations
  mat out = zeros<mat>( n_pds, n_sim ) ;
      // The matrix of outputs
  vec init( n_sim, fill::randn ) ;
  init = x0 + init * sig_eps ;
      // The initial draw is from N(0,eps)
  for( int i = 0; i < n_sim ; i++ ){
    out.col(i) = ar1_sim( n_pds, rho, sig_eps, true, init(i) ) ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat multi_norm_thresh( int n_sim, int n_pds, double rho, double sig_eps ){
// Code to create multiple random thresholds
  mat out = zeros<mat>( n_pds, n_sim ) ;
      // The matrix of outputs
  for( int i = 0; i < n_sim ; i++ ){
    out.col(i) = norm_thresh( n_pds, rho, sig_eps ) ;
  }
  return out ;
}

// [[Rcpp::export]]
List multi_thresh_filter( arma::mat x, arma::mat theta_hat, arma::imat y, 
                          arma::rowvec mu_sig2_0, double sig_eps, double rho ){
// Returns a list of the mean and variances from applying the threshold filter
// to a large sample of simulations
  int n_pds = x.n_rows ;
  int n_sim = x.n_cols ;
      // The parameters of the simulation(s)
  mat temp = zeros( n_sim, 2 ) ;
      // Temporary container for the processed signal
  mat mu = zeros( n_pds + 1, n_sim ) ;
  mat sig2 = zeros( n_pds + 1, n_sim ) ;
      // Initialize the output containers
  mat theta_hat_T = trans(theta_hat) ;
  imat y_T = trans(y) ;
      // Transpose y
  for( int i=0 ; i < n_sim ; i++ ){
    temp = thresh_filter( mu_sig2_0, theta_hat_T.row(i), sig_eps, rho, y_T.row(i) ) ;
//        Rcout << "temp:\n" << temp << std::endl ;
    mu.col(i) = temp.col(0) ;
    sig2.col(i) = temp.col(1) ;
  }
  mat mu_t = trans(mu) ;
  mat sig2_t = trans(sig2) ;
      // Format the output matrices
  List out ;
  out["mu"] = mu_t ;
  out["sig2"] = sig2_t ;
      // The output list
  return out ;
}

// [[Rcpp::export]]
List multi_gauss_filter( arma::mat x, arma::mat theta_hat, arma::imat y, 
                          arma::rowvec mu_sig2_0, double sig_eps, double rho ){
// Returns a list of the mean and variances from applying the threshold filter
// to a large sample of simulations
  int n_pds = x.n_rows ;
  int n_sim = x.n_cols ;
      // The parameters of the simulation(s)
  mat temp = zeros( n_sim, 2 ) ;
      // Temporary container for the processed signal
  mat mu = zeros( n_pds + 1, n_sim ) ;
  mat sig2 = zeros( n_pds + 1, n_sim ) ;
      // Initialize the output containers
  mat theta_hat_T = trans(theta_hat) ;
  imat y_T = trans(y) ;
      // Transpose y
  for( int i=0 ; i < n_sim ; i++ ){
    temp = gauss_filter( mu_sig2_0, theta_hat_T.row(i), sig_eps, rho, y_T.row(i) ) ;
//        Rcout << "temp:\n" << temp << std::endl ;
    mu.col(i) = temp.col(0) ;
    sig2.col(i) = temp.col(1) ;
  }
  mat mu_t = trans(mu) ;
  mat sig2_t = trans(sig2) ;
      // Format the output matrices
  List out ;
  out["mu"] = mu_t ;
  out["sig2"] = sig2_t ;
      // The output list
  return out ;
}


// To add:
// 1) mu, sigma updating
// 2) measurement of RMSE and bias @ each step