####################################################################
## charts.R
## Script to create charts for the nonlinear filter problem
## Philip Barrett, Chicago 09mar2016
####################################################################

rm(list=ls())
Rcpp::sourceCpp('momErr.cpp')
library(filters)
library(scales)
library(MASS)
library(nleqslv)

## The plot of mu and sigma updated
XX <- seq( -4, 4, length.out = 401 )
rho <- .9
sig.eps <- sqrt( 1 - rho ^ 2 )
ff <- sapply( XX, mu_sig2_update, mu_sig2=c(0,sig.eps/sqrt(1-rho^2)), sig_eps=sig.eps, rho=rho, y=1 )
gg <- sapply( XX, mu_sig2_update, mu_sig2=c(0,sig.eps/sqrt(1-rho^2)), sig_eps=sig.eps, rho=rho, y=0 )

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/mu_prime.pdf')
  plot( XX, ff[1,], ylim=c(-4,4), type='l', lwd=2, col='blue', xlab=expression(psi), 
        ylab=expression(paste(mu, "'" ) ) )
  lines( XX, gg[1,], ylim=c(-4,4), type='l', lwd=2, col='red' )
  abline( h=0, lwd=.5 )
  legend('topleft', c('y=1', 'y=0'), lwd=2, col=c('blue', 'red'), bty='n' )
dev.off()

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/sigma_prime.pdf')
  plot( XX, sqrt( ff[2,] ), type='l', lwd=2, col='blue', xlab=expression(psi), 
        ylab=expression(paste(sigma^2, "'" ) ) )
  lines( XX, sqrt( gg[2,] ), type='l', lwd=2, col='red' )
  legend('right', c('y=1', 'y=0'), lwd=2, col=c('blue', 'red'), bty='n' )
dev.off()
### Create the simulations for comparing the threshold and other filters ###
set.seed(654)
theta.hat <- 0

# The UKF parameters #
Q <- sig.eps^2
R <- .0 ^ 2
f <- function(x) rho * x
g <- function(x) if( x > theta.hat ) 1 else 0

# Create the simulation #
x.0 <- rnorm( 1, 0, sig.eps)
K <- 200
# The initial point and the length of the simulation
v.x <- c( ar1_sim( K, rho, sig.eps ) )
v.y <- 1 * ( v.x > theta.hat )

# v.y <- v.x <- rep(0, K)
# v.x[1] <- x.0
# v.y[1] <- g( v.x[1] ) + rnorm( 1, 0, R )
# for( i in 2:K ){
#   v.x[i] <- f( v.x[i-1] ) + rnorm( 1, 0, sqrt( Q ) )
#   v.y[i] <- g( v.x[i] ) + rnorm( 1, 0, sqrt( R ) )
# }   # Create the simulated state and signal

kappa <- 10
mu.sig2.0 <- c( 0, 1 )
thresh.ukf <- ukf.compute( mu.sig2.0[1], mu.sig2.0[2] , v.y, f, g, Q, R, 1, alpha=1, kappa=kappa, quad = F )
# thresh.ukf.mc <- ukf.compute( mu.sig2.0[1], mu.sig2.0[2] , v.y, f, g, Q, R, 1, alpha=1, kappa=kappa, quad = F, n.mc=10000 )
thresh.ukf.quad <- ukf.compute( mu.sig2.0[1], mu.sig2.0[2] , v.y, f, g, Q, R, 1, alpha=1, kappa=kappa, quad = T )
# The UKF (using various integration rules)
thresh <- thresh_filter( mu.sig2.0, rep(0,K), sig.eps, rho, v.y )
thresh.gf <- gauss_filter( mu.sig2.0, rep(0,K), sig.eps, rho, v.y )
# The threshold filter

#### THIS CHART INCLUDED ####
pdf('/home/philip/Dropbox//2016/Research/thesis/charts/dyn_thresh.pdf')
  plot( c(1,K), range( c( v.x, thresh[,1] + sqrt(thresh[,2]), 
                  thresh[,1] - sqrt(thresh[,2]) ) ), type='n', xlab='Period', 
        ylab='x' )
  points( 1:K, 1.02 * v.y - .01, pch=19, col=alpha('darkgreen', .5), cex=.5 )
  lines( 1:K, thresh[-(K+1),1], col='blue', lwd=2 )
  lines( 1:K, thresh[-(K+1),1] + sqrt(thresh[-(K+1),2]), col='blue', lty=2 )
  lines( 1:K, thresh[-(K+1),1] - sqrt(thresh[-(K+1),2]), col='blue', lty=2 )
  lines( 1:K, v.x, lwd=2 )
  legend( 'topright', c( 'x', 'Threshold filter mean', 
                            'Plus/minus one std dev', 'Signal' ), 
          lwd=c(2,2,1,0), lty=c(1,1,2, NA), pch=c(NA,NA,NA,19), bty='n',
          col=c( 'black','blue', 'blue', alpha( 'darkgreen', .5) ))
  abline( h=0, lwd=.5 )
dev.off()

mu.sig.bar.fun <- function( mu.sig.bar, y ){
  out <- mu.sig.bar - mu_sig2_update( mu.sig.bar, theta.hat, sig.eps, rho, y )
}

mu.sig.bar.1 <- nleqslv( c( 1, .5 ), mu.sig.bar.fun, y=1 )
mu.sig.bar.0 <- nleqslv( c( 1, .5 ), mu.sig.bar.fun, y=0 )

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/xsect_thresh.pdf')
  plot( thresh[-(1:5),1], thresh[-(1:5),2], xlab=expression(mu), ylab=expression(sigma^2),
        pch=19, col='blue', cex=.5, xlim=c(-1,1), ylim=c(.2, .5) )
  points( c( mu.sig.bar.1$x[1], mu.sig.bar.0$x[1] ), 
          c( mu.sig.bar.1$x[2], mu.sig.bar.0$x[2] ), pch=19 )
  legend( 'bottomright', c('Mean-variance pairs', 'Limit point'), pch=19,
          col=c('blue','black'), bty='n' )
dev.off()

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/dyn_gf.pdf')
  plot( c(1,K), range( c( v.x, thresh.gf[,1] + sqrt(thresh.gf[,2]), 
                          thresh.gf[,1] - sqrt(thresh.gf[,2]) ) ), type='n', xlab='Period', 
        ylab='x' )
  lines( 1:K, thresh[-(K+1),1], col='blue', lwd=2 )
  points( 1:K, 1.02 * v.y - .01, pch=19, col=alpha('darkgreen', .5), cex=.5 )
  lines( 1:K, thresh.gf[-(K+1),1], col='red', lwd=2 )
  lines( 1:K, thresh.gf[-(K+1),1] + sqrt(thresh.gf[-(K+1),2]), col='red', lty=2 )
  lines( 1:K, thresh.gf[-(K+1),1] - sqrt(thresh.gf[-(K+1),2]), col='red', lty=2 )
  lines( 1:K, v.x, lwd=2 )
  legend( 'topright', c( 'x', 'Exact Gaussian filter mean', 
                         'Plus/minus one std dev', 'Threshold filter mean', 'Signal' ), 
          lwd=c(2,2,1,2,0), lty=c(1,1,2,1, NA), pch=c(NA,NA,NA,NA,19), bty='n',
          col=c( 'black','red', 'red', 'blue', alpha( 'darkgreen', .5) ))
  abline( h=0, lwd=.5 )
dev.off()

pdf('/home/philip/Dropbox//2016/Research/thesis/charts/xsect_gf.pdf')
  plot( thresh.gf[-(1:20),1], thresh.gf[-(1:20),2], xlab=expression(mu), ylab=expression(sigma^2),
        pch=19, col='red', cex=.5, xlim=c(-1,1), ylim=c(.2, .5) )
#   points( c( mu.sig.bar.1$x[1], mu.sig.bar.0$x[1] ), 
#           c( mu.sig.bar.1$x[2], mu.sig.bar.0$x[2] ), pch=19 )
  legend( 'bottomright', c('Mean-variance pairs', 'Limit point'), pch=19,
          col=c('red','black'), bty='n' )
dev.off()


plot( c(1,K), range( c(thresh.ukf$m, v.x, thresh[,1]) ), type='n', xlab='Period', 
      ylab='x' )
points( 1:K, 1.1 * sd(v.x) * ( 2*v.y-1 ), pch=19, col=alpha('darkgreen', .5), cex=.5 )
lines( 1:K, thresh.gf[-(K+1),1], col='red', lwd=2 )
lines( 1:K, thresh.ukf$m.pred[-(K+1)], col='red', lwd=1, lty=2 )
# lines( 1:K, thresh.ukf.mc$m.pred[-(K+1)], col='red', lwd=1, lty=3 )
# lines( 1:K, thresh.ukf.quad$m.pred[-(K+1)], col='red', lwd=1, lty=3 )
# First point is the period 0 predictor for period 1 => Last point predicts
# period K+1
lines( 1:K, thresh[-(K+1),1], col='blue', lwd=2 )
# Likewise
# lines( 1:K, thresh.ukf$m + sqrt( c( thresh.ukf$P.pred[-K] ) ), col='red', lwd=2, lty=2 )
# lines( 1:K, thresh.ukf$m - sqrt( c( thresh.ukf$P.pred[-K] ) ), col='red', lwd=2, lty=2 )
lines( 1:K, v.x, lwd=2 )
legend( 'bottomright', c( 'x', 'Threshold filter', 'Exact Gaussian Filter',
                          'Unscented Kalman Filter', 'Signal' ), 
        lwd=c(2,2,2,1,1,0), lty=c(1,1,1,2,3, NA), pch=c(NA,NA,NA,NA,NA,19), bty='n',
        col=c( 'black','blue', 'red', 'red', 'red', alpha( 'darkgreen', .5) ))
abline( h=0, lwd=.5 )

plot( 1:K, sqrt(thresh[-(K+1),2]), type='l', lwd=2, col='red' )
lines( 1:K, sqrt(thresh[-(K+1),2]), type='l', lwd=2, col='blue' )

rmse <- sqrt( cumsum( ( v.x - thresh[-(K+1),1] ) ^ 2 ) / 1:K )
rmse.gf <- sqrt( cumsum( ( v.x - thresh.gf[-(K+1),1] ) ^ 2 ) / 1:K )
rmse.ukf <- sqrt( cumsum( ( v.x - thresh.ukf$m.pred[-(K+1)] ) ^ 2 ) / 1:K )
plot( c(1,K), range( rmse, rmse.gf ), type='n', xlab='Period', ylab='Rolling RMSE' )
lines( 1:K, rmse, col='blue', lwd=2)
lines( 1:K, rmse.gf, col='red', lwd=2 )

bias <- cumsum( ( v.x - thresh[-(K+1),1] ) ) / 1:K
bias.gf <- cumsum( ( v.x - thresh.gf[-(K+1),1] ) ) / 1:K
bias.ukf <- cumsum( v.x - thresh.ukf$m.pred[-(K+1)] ) / 1:K 
plot( c(1,K), range( bias, bias.gf ), type='n', xlab='Period', ylab='Rolling bias' )
lines( 1:K, bias, col='blue', lwd=2)
lines( 1:K, bias.gf, col='red', lwd=2 )
abline( h=0, lwd=.5 )

#### Now generate a bunch of simulations and see the properties of the errors ###
set.seed(4321)
n.sim <- 100000
n.pds <- 20
multi.x <- multi_ar1_sim( n.sim, n.pds, rho, 0, sig.eps )
multi.theta.hat <- 0.0 * multi.x # multi_norm_thresh( n.sim, n.pds, rho, sig.eps )
multi.y <- 1 * ( multi.x > multi.theta.hat )
multi.thresh <- multi_thresh_filter( multi.x, multi.theta.hat, multi.y, 
                                     c( 0, sig.eps^2 ), sig.eps, rho )
multi.gauss <- multi_gauss_filter( multi.x, multi.theta.hat, multi.y, 
                                     c( 0, sig.eps^2 ), sig.eps, rho )
err <- multi.thresh$mu[,-(n.pds+1)] - t( multi.x )
bias <- apply( err, 2, mean )
rmse <- apply( err, 2, sd )
mse <- apply( err, 2, var )
err.gf <- multi.gauss$mu[,-(n.pds+1)] - t( multi.x )
bias.gf <- apply( err.gf, 2, mean )
rmse.gf <- apply( err.gf, 2, sd )
mse.gf <- apply( err.gf, 2, var )

sig.mean <- apply( sqrt( multi.thresh$sig2 ), 2, mean )

# multi.thresh.ukf <- list( m.pred=0*multi.thresh$mu, P.pred=0*multi.thresh$sig2 )
# for( i in 1:n.sim ){
#   temp <- ukf.compute( 0, sig.eps^2, multi.y[,i], f, g, Q, R, 1, 
#                        alpha=1, kappa=kappa, quad = F )
#   multi.thresh.ukf$m.pred[i,] <- temp$m.pred
#   multi.thresh.ukf$P.pred[i,] <- temp$P.pred
# }
# err.ukf <- multi.thresh.ukf$m.pred[,-(n.pds+1)] - t( multi.x )
# bias.ukf <- apply( err.ukf, 2, mean )
# rmse.ukf <- apply( err.ukf, 2, sd )
# mse.ukf <- apply( err.gf, 2, var )

#### THIS CHART INCLUDED ####
# plot( 1:n.pds, rmse.ukf, col='red', lty=2, lwd=2, type='l', xlab='Periods', ylab='RMSE' )
plot( 1:n.pds, rmse.gf, col='red', lwd=2, type='l', xlab='Periods', ylab='RMSE' )
# lines( 1:n.pds, rmse.gf, col='red', lwd=2 )
lines( 1:n.pds, rmse, col='blue', lwd=2 )
# lines( 1:20, sqrt(apply(multi.gauss$sig2[,-(n.pds+1)],2,mean)), lty=2, col='red' )

plot( 1:20, apply(multi.x, 1, sd), lwd=2, type='l', xlab='Periods', 
      ylab='State sd' )
tot.var.thresh <- apply(multi.thresh$sig2[,-(n.pds+1)],2,mean) + apply(multi.thresh$mu[,-(n.pds+1)],2,var)
tot.var.gf <- apply(multi.gauss$sig2[,-(n.pds+1)],2,mean) + apply(multi.gauss$mu[,-(n.pds+1)],2,var)
lines( 1:20, sqrt(tot.var.thresh), lwd=2, col='blue' )
lines( 1:20, sqrt(tot.var.gf), lwd=2, col='red' )
legend( 'bottomright', c('State variance', 'Total variance: Threshold filter', 'Total variance: Exact Gaussian filter'),
        bty='n', lwd=2, col=c( 'black', 'blue', 'red' ) )

plot( 1:20, 1 - tot.var.gf / apply(multi.x, 1, var), lwd=2, col='red', type='l' )
lines( 1:20, 1 - tot.var.thresh / apply(multi.x, 1, var), lwd=2, col='blue' )


#### NOW DO CONDITIONAL BIAS CHARTS ####
n.pds <- 100000
burn <- 1000
x.lr <- c( ar1_sim( n.pds + burn, rho, sig.eps ) )
    # Long run x
theta.hat.lr <- 0 * x.lr
y.lr <- 1 * ( x.lr > theta.hat.lr )

# Create the filters
thresh.lr <- thresh_filter( c(0,sig.eps^2), theta.hat.lr, sig.eps, rho, y.lr )
thresh.lr.gf <- gauss_filter( c(0,sig.eps^2), theta.hat.lr, sig.eps, rho, y.lr )
thresh.lr.ukf <- ukf.compute( 0, sig.eps^2, y.lr, f, g, Q, R, 1, alpha=1, kappa=kappa )

# De-burn
thresh.lr <- thresh.lr[-(1:burn),]
thresh.lr.gf <- thresh.lr.gf[-(1:burn),]
m.thresh.lr.ukf <- cbind( thresh.lr.ukf$m.pred[-(1:burn)], 
                          thresh.lr.ukf$P.pred[-(1:burn)] )
x.lr <- x.lr[-(1:burn)]
y.lr <- y.lr[-(1:burn)]

# Create the conditional biases
bias.pos.lr <- mean( thresh.lr[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1] )
bias.pos.lr.gf <- mean( thresh.gf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1] )
bias.pos.lr.ukf <- mean( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1] )
bias.neg.lr <- mean( thresh.lr[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0] )
bias.neg.lr.gf <- mean( thresh.lr.gf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0] )
bias.neg.lr.ukf <- mean( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0] )
n.same <- sequence(rle(y.lr)$lengths)
    # The number of identical signals
table( n.same, y.lr )
bias.p.seq <- c( by( thresh.lr[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                  n.same[y.lr==1], mean ) )[-1]
bias.p.seq.gf <- c( by( thresh.lr.gf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                     n.same[y.lr==1], mean ) )[-1]
bias.p.seq.ukf <- c( by( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                  n.same[y.lr==1], mean ) )[-1]
bias.n.seq <- c( by( thresh.lr[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                  n.same[y.lr==0], mean ) )[-1]
bias.n.seq.gf <- c( by( thresh.lr.gf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                     n.same[y.lr==0], mean ) )[-1]
bias.n.seq.ukf <- c( by( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                  n.same[y.lr==0], mean ) )[-1]

rmse.p.seq <- c( by( thresh.lr[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                     n.same[y.lr==1], function(x) sqrt(mean(x^2)) ) )[-1]
rmse.p.seq.gf <- c( by( thresh.lr.gf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                     n.same[y.lr==1], function(x) sqrt(mean(x^2)) ) )[-1]
rmse.p.seq.ukf <- c( by( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==1] - x.lr[y.lr==1], 
                         n.same[y.lr==1], function(x) sqrt(mean(x^2)) ) )[-1]
rmse.n.seq <- c( by( thresh.lr[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                     n.same[y.lr==0], function(x) sqrt(mean(x^2)) ) )[-1]
rmse.n.seq.gf <- c( by( thresh.lr.gf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                     n.same[y.lr==0], function(x) sqrt(mean(x^2)) ) )[-1]
rmse.n.seq.ukf <- c( by( m.thresh.lr.ukf[-(n.pds+1),1][y.lr==0] - x.lr[y.lr==0], 
                         n.same[y.lr==0], function(x) sqrt(mean(x^2)) ) )[-1]


## Now plot them
plot( c(1,10), range( bias.p.seq[1:10], bias.p.seq.ukf[1:10], 
                      bias.n.seq[1:10], bias.n.seq.ukf[1:10] ), type='n' )
lines( 1:10, bias.p.seq[1:10], lwd=2, col='blue' )
lines( 1:10, bias.n.seq[1:10], lwd=2, col='blue', lty=2 )
lines( 1:10, bias.p.seq.gf[1:10], lwd=2, col='red' )
lines( 1:10, bias.n.seq.gf[1:10], lwd=2, col='red', lty=2 )
lines( 1:10, bias.p.seq.ukf[1:10], col='red' )
lines( 1:10, bias.n.seq.ukf[1:10], col='red', lty=2 )
abline(h=0, lwd=.5)

plot( c(1,10), range( 0, rmse.p.seq[1:10], rmse.p.seq.ukf[1:10], 
                      rmse.n.seq[1:10], rmse.n.seq.ukf[1:10] ), type='n' )
lines( 1:10, rmse.p.seq[1:10], lwd=2, col='blue' )
# lines( 1:10, rmse.n.seq[1:10], lwd=2, col='blue', lty=2 )
lines( 1:10, rmse.p.seq.ukf[1:10], lwd=2, col='red' )
# lines( 1:10, rmse.n.seq.ukf[1:10], lwd=2, col='red', lty=2 )
abline(h=0, lwd=.5)
