


rm(list=ls())
Rcpp::sourceCpp('momErr.cpp')

exact.1.pd.cdf.fn <- function( v.x, mu.sig2, x.hat, y, rho, sig.eps ) {
  mu <- mu.sig2[1]
  sig2 <- mu.sig2[2]
  sig <- sqrt(sig2)
  gamma2 <- rho ^ 2 * sig2 + sig.eps
  gamma <- sqrt( gamma2 )
      # Extract
  A <- pnorm( rho * sig2 * v.x, (gamma2 * x.hat - mu * sig.eps ^ 2 ) , ( sig * sig.eps * gamma ) )
  B <- pnorm( mu, x.hat, sig )
  if( y == 0 ){
    A <- 1 - A
    B <- 1 - B
  } 
  C <- dnorm( v.x, rho * mu, gamma )
      # Components of the pdf
  pdf <- A / B * C
  d <- diff( v.x )
  cdf <- cumsum(pdf[-1]*d)
      # The cdf
  return( cdf / tail(cdf,1) )
}

exact.err.fn <- function( rho, n.pds, n.markov, burn ){
  message('rho = ', rho)
  sig.eps <- sqrt( 1 - rho ^ 2 )
  sim <- ar1_sim( n.pds, rho, sig.eps )
      # The simulation
  x.hat <- norm_thresh( n.pds, rho, sig.eps )
  y <- signal_create( sim, x.hat )
      # The signals
  v.x <- seq( -6, 6, length.out = 501 )
  thresh <- thresh_filter( c(0,1), x.hat, sig.eps, rho, y )
  gf <- gauss_filter( c(0,1), x.hat, sig.eps, rho, y )
      # The two filters
  thresh.cdf <- t( apply( thresh, 1, function(ms) pnorm( v.x, ms[1], sqrt(ms[2]) ) ) )
  gf.cdf <- t( apply( gf, 1, function(ms) pnorm( v.x, ms[1], sqrt(ms[2]) ) ) )
      # The filter cdfs
  exact.1.pd.cdf <- t( sapply( 1:(nrow(thresh)-1), 
                               function(i) exact.1.pd.cdf.fn( v.x, thresh[i,], x.hat[i], y[i], rho, sig.eps ) ) )
      # The exact cdf
  exact.1.pd.pdf <- t( apply( exact.1.pd.cdf, 1, function(x) c( x[1], diff(x) ) ) ) # / diff(v.x) ) )
      # The pdf (use for integration)
  thresh.cdf.err <- apply( exact.1.pd.pdf * abs( thresh.cdf[-1,-1] - exact.1.pd.cdf ), 1, sum )
  gf.cdf.err <- apply( exact.1.pd.pdf * abs( gf.cdf[-1,-1] - exact.1.pd.cdf ), 1, sum )
      # The errors
  cdf.tab <- t( apply( cbind( threshold=thresh.cdf.err, gaussian=gf.cdf.err ), 2, 
                       function( x ) c( mean=mean(x), max=max(x) ) ) )
  return( list( rho=rho, cdf=cdf.tab ) )
}

n.pds <- 20100
n.markov <- 600
burn <- 1000
n.rho <- 40

v.rho=seq( 0, .975, length.out=n.rho )

l.err <- lapply( v.rho, exact.err.fn, n.pds=n.pds, n.markov=n.markov, burn=burn )

save( l.err, file='err_exact.rdata' )

max.cdf <- sapply( l.err, function(x) c( x$rho, x$cdf[,'mean'] ) )
mean.cdf <- sapply( l.err, function(x) c( x$rho, x$cdf[,'mean'] ) )

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/exact_cdf_err.pdf')
plot( max.cdf[1,], max.cdf['gaussian',], type='l', xlim=c(0,1), 
      xlab=expression(rho), ylab='CDF error', col='red', lwd=2 )
lines( max.cdf[1,], max.cdf['threshold',], type='l', col='blue', lwd=2 )
lines( max.cdf[1,], mean.cdf['threshold',], type='l', col='blue', lwd=2, lty=2 )
lines( max.cdf[1,], mean.cdf['gaussian',], type='l', col='red', lwd=2, lty=2 )
legend( 'topleft', c('Threshold filter max error', 'Threshold filter mean error', 
                     'Exact Gaussian filter max error', 'Exact Gaussian filter mean error' ), 
        lwd=2, lty=c(1,2,1,2), col=c('blue', 'blue', 'red', 'red'), bty='n' )
dev.off()


