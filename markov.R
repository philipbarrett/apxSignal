

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
  thresh.err <- apply( abs( cbind( thresh[-1,], 0 ) - markov.mom )[-(1:burn),], 2, mean ) / markov.mom.ave
  gf.err <- apply( abs( cbind( gf[-1,], 0 ) - markov.mom )[-(1:burn),], 2, mean ) / markov.mom.ave
      # Error on the two filters relative to the distribution sd
  
  markov.cdf <- t(apply( markov$filter, 1, cumsum ))[ -(1:burn), ]
      # The Markov cdf
  thresh.cdf <- t( apply( thresh[-1,], 1, function(x) 
    pnorm( markov$b[-1], x[1], sqrt( x[2] ) ) ) )[ -(1:burn), ]
  gf.cdf <- t( apply( gf[-1,], 1, function(x) 
    pnorm( markov$b[-1], x[1], sqrt( x[2] ) ) ) )[ -(1:burn), ]
      # Compute the cdfs
  thresh.cdf.err <- apply( markov$filter[-(1:burn),] * abs( markov.cdf - thresh.cdf ), 1, mean )
  gf.cdf.err <- apply( markov$filter[-(1:burn),] * abs( markov.cdf - gf.cdf ), 1, mean )
    # Compute the errors
  
  cdf.tab <- t( apply( cbind( threshold=thresh.cdf.err, gaussian=gf.cdf.err ), 2, 
                       function( x ) c( mean=mean(x), max=max(x) ) * 100 ) )
  mom.tab <- rbind( t(thresh.err), t(gf.err) )
  rownames( mom.tab) <- rownames( cdf.tab )
      # Create the error tables
  return( list( rho=rho, mom=mom.tab, cdf=cdf.tab ) )
}

mu <- 0
rho <- .95
sig.eps <- sqrt( 1 - rho ^ 2 )
n.pds <- 10100
n.markov <- 500
burn <- 1000

v.rho=seq( 0, .95, length.out=20 )
l.err <- mclapply( v.rho, this.err.markov, n.pds=n.pds, n.markov=n.markov, burn=burn )

max.cdf <- sapply( l.err, function(x) c( x$rho, x$cdf[,'max'] ) )
plot( max.cdf[1,], max.cdf['gaussian',], type='l' )
lines( max.cdf[1,], max.cdf['threshold',], type='l', col=2 )
