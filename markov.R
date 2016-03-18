
rm(list=ls())
Rcpp::sourceCpp('momErr.cpp')
library(parallel)

this.err.markov <- function( rho, n.pds, n.markov, burn ){
  message('rho = ', rho)
  sig.eps <- sqrt( 1 - rho ^ 2 )
  markov <- markov_filter_sim( mu, sig.eps, rho, n.markov, n.pds )
  thresh <- thresh_filter( c( 0, 1 ), markov$x_hat, sig.eps, rho, markov$y )
  gf <- gauss_filter( c( 0, 1 ), markov$x_hat, sig.eps, rho, markov$y )
  markov.mu <- markov$filter %*% markov$m
  markov.mom <- cbind( mu=c(markov.mu), 
                       sig2=apply( markov$filter * ( rep(1,n.pds) %*% t(markov$m) - 
                                                       markov.mu %*% t(rep(1,n.markov)) ) ^ 2, 1, sum ), 
                       skew3=( apply( markov$filter * ( rep(1,n.pds) %*% t(markov$m) - 
                                                          markov.mu %*% t(rep(1,n.markov)) ) ^ 3, 1, sum ) ) )
      # Create the moments of the various filters
  
  markov.sig.ave <- sapply( mean( markov.mom[,'sig2'] ), function(x) c( mu=x^.5, sig2=x, skew=x ^ (3/2) ) )
  thresh.err <- apply( abs( cbind( thresh[-1,], 0 ) - markov.mom )[-(1:burn),], 2, mean ) / markov.sig.ave
  gf.err <- apply( abs( cbind( gf[-1,], 0 ) - markov.mom )[-(1:burn),], 2, mean ) / markov.sig.ave
      # Error on the two filters relative to the distribution sd
  
  markov.cdf <- t(apply( markov$filter, 1, cumsum ))[ -(1:burn), ]
      # The Markov cdf
  thresh.cdf <- t( apply( thresh[-1,], 1, function(x) 
    pnorm( markov$b[-1], x[1], sqrt( x[2] ) ) ) )[ -(1:burn), ]
  gf.cdf <- t( apply( gf[-1,], 1, function(x) 
    pnorm( markov$b[-1], x[1], sqrt( x[2] ) ) ) )[ -(1:burn), ]
      # Compute the cdfs
  thresh.cdf.err <- apply( markov$filter[-(1:burn),] * abs( markov.cdf - thresh.cdf ), 1, sum )
  gf.cdf.err <- apply( markov$filter[-(1:burn),] * abs( markov.cdf - gf.cdf ), 1, sum )
    # Compute the errors
  
  cdf.tab <- t( apply( cbind( threshold=thresh.cdf.err, gaussian=gf.cdf.err ), 2, 
                       function( x ) c( mean=mean(x), max=max(x) ) ) )
  mom.tab <- rbind( t(thresh.err), t(gf.err) )
  rownames( mom.tab) <- rownames( cdf.tab )
      # Create the error tables
  return( list( rho=rho, mom=mom.tab, cdf=cdf.tab ) )
}

mu <- 0
n.pds <- 20100
n.markov <- 600
burn <- 1000
n.rho <- 40

v.rho=seq( 0, .975, length.out=n.rho )
l.err <- lapply( v.rho, this.err.markov, n.pds=n.pds, n.markov=n.markov, burn=burn )
# l.err <- mclapply( v.rho, this.err.markov, n.pds=n.pds, n.markov=n.markov, burn=burn )

save( l.err, file='err_markov.rdata' )

max.cdf <- sapply( l.err, function(x) c( x$rho, x$cdf[,'max'] ) )
mean.cdf <- sapply( l.err, function(x) c( x$rho, x$cdf[,'mean'] ) )

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/markov_cdf_err.pdf')
plot( max.cdf[1,], max.cdf['gaussian',], type='l', xlim=c(0,1), 
      xlab=expression(rho), ylab='CDF error', col='red', lwd=2 )
lines( max.cdf[1,], max.cdf['threshold',], type='l', col='blue', lwd=2 )
lines( max.cdf[1,], mean.cdf['threshold',], type='l', col='blue', lwd=2, lty=2 )
lines( max.cdf[1,], mean.cdf['gaussian',], type='l', col='red', lwd=2, lty=2 )
legend( 'topleft', c('Threshold filter max error', 'Threshold filter mean error', 
                     'Exact Gaussian filter max error', 'Exact Gaussian filter mean error' ), 
        lwd=2, lty=c(1,2,1,2), col=c('blue', 'blue', 'red', 'red'), bty='n' )
dev.off()

mu.err <- sapply( l.err, function(x) c( x$rho, x$mom[,'mu'] ) )
pdf('/home/philip/Dropbox//2016/Research/thesis/charts/markov_mu_err.pdf')
plot( mu.err[1,], mu.err['gaussian',], type='l', xlim=c(0,1), 
      xlab=expression(rho), ylab=expression( paste( 'Ave abs error: ', mu ) ), col='red', lwd=2 )
lines( mu.err[1,], mu.err['threshold',], type='l', col='blue', lwd=2 )
legend( 'topleft', c('Threshold filter', 'Exact Gaussian filter'), lwd=2, 
        col=c('blue', 'red'), bty='n' )
dev.off()

sig.err <- sapply( l.err, function(x) c( x$rho, sqrt( x$mom[,'sig2'] ) ) )
pdf('/home/philip/Dropbox//2016/Research/thesis/charts/markov_mu_err.pdf')
plot( sig.err[1,], sig.err['gaussian',], type='l', xlim=c(0,1), 
      xlab=expression(rho),  ylab=expression( paste( 'Ave abs error: ', sigma ) ), col='red', lwd=2 )
lines( sig.err[1,], sig.err['threshold',], type='l', col='blue', lwd=2 )
legend( 'topleft', c('Threshold filter', 'Exact Gaussian filter'), lwd=2, 
        col=c('blue', 'red'), bty='n' )
dev.off()
