######################################################################
## Bayesian recursive smoothing via EpiFilter
# From: Parag, KV (2021). Improved real-time estimation of reproduction numbers at
# low case incidence and between epidemic waves. PLOS Comput Biol 17(9):e1009347.
######################################################################

# Assumptions
# - observation model is Poisson growth rate equation
# - reproduction number state space model is a simple diffusion
# - must have run epiFilter first to obtain forward distribution pr
# - method makes a backward pass to generate qr

# Inputs - grid on reproduction numbers (rgrid), size of grid (m), filtered posterior (pr),
# update pre-filter (prup), max time (nday), state transition matrix (pstate), confidence (a)

# Output - mean (rmean), median (rmed), lower (rlow) amd upper (rhigh) quantiles of estimates,
# smoothed posterior over r (qr) which is backwards and forwards

epiSmootherGrowth <- function(rgrid, m, pr, prup, nday, pstate, a){
  
  # Last smoothed distribution same as filtered
  qr = matrix(0, nday, m); qr[nday, ] = pr[nday, ]
  
  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # remove zeros
    prup[i+1, prup[i+1, ] == 0] = 10^-8
    
    # Integral term in smoother
    integ = qr[i+1, ]/prup[i+1, ]
    integ = integ%*%pstate
    
    # Smoothed posterior over rgrid
    qr[i, ] = pr[i, ]*integ
    # Force a normalisation
    qr[i, ] = qr[i, ]/sum(qr[i, ]);
  }
  
  # Mean, median estimates of r and variance
  rmean = rep(0, nday); rmed = rmean; rvar = rmean
  # 50% and 95% (depends on a) confidence on r
  rci = matrix(0, 4, nday)
  
  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean, var and CDF
    rmean[i] = qr[i, ]%*%rgrid
    rvar[i] = qr[i, ]%*%(rgrid^2) - rmean[i]^2
    rcdf = cumsum(qr[i, ])
    
    # Quantiles for estimates
    idm = which(rcdf >= 0.5); rmed[i] = rgrid[idm[1]]
    id1 = which(rcdf >= a, 1); id2 = which(rcdf >= 1-a, 1)
    id3 = which(rcdf >= 0.25, 1); id4 = which(rcdf >= 0.75, 1)
    rci[1, i] = rgrid[id1[1]]; rci[2, i] = rgrid[id2[1]]
    rci[3, i] = rgrid[id3[1]]; rci[4, i] = rgrid[id4[1]]
  }
  
  # Main outputs: estimates of r and states
  out = list(rmed, rci, rmean, qr, rvar)
  names(out) = c('rmed', 'rci', 'rmean', 'qr', 'rvar')
  epiSmootherGrowth = out
}