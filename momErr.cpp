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

// [[Rcpp::export]]
double trunc_mean( double a, double b, double mu, double sig ){
// Computes the mean of a normal truncated to lie in [a,b]
  boost::math::normal norm_dist( 0, 1 ) ;
      // A normal distribution
  double a_norm = ( a - mu ) / sig ;
  double b_norm = ( b - mu ) / sig ;
      // The normalized limits
  double out = mu + sig * ( pdf( norm_dist, a_norm ) - pdf( norm_dist, b_norm ) ) / 
                              ( cdf( norm_dist, b_norm ) - cdf( norm_dist, a_norm ) ) ;
  return out ;
}

// [[Rcpp::export]]
List norm_split( double mu, double sig, int n, bool spread=false ){
// Creates the splits for a discretized normal
  vec quantiles(n+1) ;
  for( int i = 0 ; i < n + 1 ; i++ ) quantiles(i) = (double)(i) / n ;
      // The quantiles
  if( spread ) quantiles = quantiles - .15 * sin( quantiles * 2 * datum::pi ) ;
      // Use sin curve to spread the quantile targets
  vec b(n+1);
  boost::math::normal norm_dist( mu, sig ) ;
  b(0) = -datum::inf ;
  b(n) = datum::inf ;
  for( int i = 1 ; i < n ; i++ ) b(i) = quantile( norm_dist, quantiles(i) ) ;
      // The quantile values
  vec m(n) ;
  for( int i = 0 ; i < n ; i++ ) m(i) = trunc_mean( b(i), b(i+1), mu, sig ) ;
      // The conditional means
  List out ;
  out["q"] = quantiles ;
  out["m"] = m ;
  out["b"] = b ;
  return out ;
}

// [[Rcpp::export]]
arma::rowvec pt_evolve( double pt, double eps, double rho, double mu, arma::vec b ){
// Evolves a single mass-point through the AR(1) when the cutoffs for the
// discretization are given by b.  Returns a vector of weights in each b range
  boost::math::normal norm_dist( rho * pt + ( 1 - rho ) * mu, eps ) ;
      // The normal distribution for the evolution of the point
  int n = b.n_elem - 1 ;
  rowvec out(n) ;
  for( int i=0 ; i < n ; i++ ) out(i) = cdf( norm_dist, b(i+1) ) - cdf( norm_dist, b(i) ) ;
      // Fill in the weights
  return out ;
}

// [[Rcpp::export]]
List markov_disc( double mu, double eps, double rho, int n ){
// Discretize an AR(1)
  double sig = eps / sqrt( 1 - pow( rho, 2 ) ) ;
      // The unconditional variance
  List out = norm_split( mu, sig, n, true ) ;
      // Initialize the output
  vec m = out["m"] ;
  vec b = out["b"] ;
      // Extra output
  mat trans(n,n) ;
  for( int i = 0 ; i < n ; i++ ) trans.row(i) = pt_evolve( m(i), eps, rho, mu, b ) ;
      // The transition matrix
  out["trans"] = trans ;
  return out ;
}

// [[Rcpp::export]]
arma::rowvec markov_filter( arma::mat trans, arma::vec b, arma::rowvec dist, double x_hat, int y ){
// Computes the evolution of the threshold problem with a high-state markov approximation
  int n = b.n_elem - 1 ;
      // The size of the approx
  rowvec trunc(n) ;
  if( y == 1 ){
    for( int i = 0 ; i < n ; i ++ ){
      if( b(i+1) < x_hat ){
        trunc(i) = 0 ;
      }else{
        if( b(i) > x_hat ){
          trunc(i) = dist(i) ;
        }else{
          trunc(i) = dist(i) * ( b(i+1) - x_hat ) / ( b(i+1) - b(i) ) ;
        }
      }
    }
  }
  if( y == 0 ){
    for( int i = 0 ; i < n ; i ++ ){
      if( b(i+1) < x_hat ){
        trunc(i) = dist(i) ;
      }else{
        if( b(i) > x_hat ){
          trunc(i) = 0 ;
        }else{
          trunc(i) = dist(i) * ( x_hat - b(i) ) / ( b(i+1) - b(i) ) ;
        }
      }
    }
  }
  trunc = trunc / sum( trunc ) ;
      // Rescale
  rowvec out = trunc * trans ;
  return out ;
      // The answer
}

// [[Rcpp::export]]
List markov_filter_sim( double mu, double sig_eps, double rho, int n_markov, int n_pds ){
// Simulates an AR(1) process and returns the Markov filter
  vec sim = ar1_sim( n_pds, rho, sig_eps ) ;
      // The simulation
  vec x_hat = norm_thresh( n_pds, rho, sig_eps ) ;
      // Generate the thrsholds
  ivec y = signal_create( sim, x_hat ) ;
      // The sequence of signals
  mat markov = zeros( n_pds, n_markov ) ;
      // Initialize the matrix of probabilities
  List markov_def = markov_disc( mu, sig_eps, rho, n_markov ) ;
      // The definition of the Markov filter
  mat trans = markov_def["trans"] ;
  vec b = markov_def["b"] ;
  vec m = markov_def["m"] ;
  vec q = markov_def["q"] ;
      // Extract elements of the markov_def object
  vec ergodic = q.tail(n_markov) - q.head(n_markov) ;
  rowvec ergodic_r = conv_to<rowvec>::from(ergodic) ;
      // Initialize with the ergodic dstribution
  markov.row(0) = ergodic_r ;
      // The first row of the Markov filter
  for( int i = 1 ; i < n_pds ; i++ ){
    markov.row(i) = markov_filter( trans, b, markov.row(i-1), x_hat(i), y(i) ) ;
  }   // Fill in the markov filter
  List out ;
  out["sim"] = sim ;
  out["x_hat"] = x_hat ;
  out["y"] = y ;
  out["b"] = b ;
  out["q"] = q ;
  out["m"] = m ;
  out["filter"] = markov ;
  out["ergodic"] = ergodic ;
      // Format the output
  return out ;
}
