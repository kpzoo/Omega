######################################################################
## Bayesian recursive filtering - intermediate step of EpiFilter for growth rates
# From: Parag, KV (2021). Improved real-time estimation of reproduction numbers at
# low case incidence and between epidemic waves. PLOS Comput Biol 17(9):e1009347.
######################################################################

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion

# Inputs - grid on growth rates (rgrid), size of grid (m), diffusion noise (eta),
# prior on r (pr0), max time (nday), incidence or some convolution (Iday), confidence (a)

# Output - mean (rmean), median (rmed), 50% and 95% quantiles of estimates (rhat),
# causal posterior over r (pr), pre-update (prup) and state transition matrix (pstate)

# Sources of errors 
# - if rgrid is poorly specified (m too small or rmin/rmax too tight) due to rcdf not being computable
# - if Iday[1] = 0 as the epidemic must be seeded with at least one case (this is a renewal model trait)

epiFilterGrowth <- function(rgrid, m, eta, pr0, nday, Iday, a){
  
  # Probability vector for r and prior
  pr = matrix(0, nday, m); prup = pr
  pr[1, ] = pr0; prup[1, ] = pr0
  
  # Mean and median estimates
  rmean = rep(0, nday); rmed = rmean
  # 50% and 95% (depends on a) confidence on r
  rci = matrix(0, 4, nday)
  
  # Initialise mean and CDF of prior 
  rmean[1] = pr[1, ]%*%rgrid; rcdf0 = cumsum(pr0)
  
  # Initialise quartiles
  idm = which(rcdf0 >= 0.5, 1); rmed[1] = rgrid[idm[1]]
  id1 = which(rcdf0 >= a, 1); id2 = which(rcdf0 >= 1-a, 1)
  id3 = which(rcdf0 >= 0.25, 1); id4 = which(rcdf0 >= 0.75, 1)
  rci[1, 1] = rgrid[id1[1]]; rci[2, 1] = rgrid[id2[1]]
  rci[3, 1] = rgrid[id3[1]]; rci[4, 1] = rgrid[id4[1]]
  
  # Precompute state distributions for r transitions
  pstate = matrix(0, m, m);
  for(j in 1:m){
    pstate[j, ] = dnorm(rgrid[j], rgrid, eta)
  }
  pstate = t(pstate)
  
  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    if(Iday[i-1] == 0){rate = exp(rgrid)}else{rate = Iday[i-1]*exp(rgrid)} # can setup as avg of days
    # Probabilities of observations
    pI = dpois(Iday[i], rate)
    
    # State predictions for r
    prup[i, ]  = pr[i-1, ]%*%pstate
    # Update to posterior over r
    pr[i, ] = prup[i, ]*pI
    pr[i, ] = pr[i, ]/sum(pr[i, ])
    
    # Posterior mean and CDF
    rmean[i] = pr[i, ]%*%rgrid
    rcdf = cumsum(pr[i, ])
    
    # Quantiles for estimates
    idm = which(rcdf >= 0.5, 1); rmed[i] = rgrid[idm[1]]
    id1 = which(rcdf >= a, 1); id2 = which(rcdf >= 1-a, 1)
    id3 = which(rcdf >= 0.25, 1); id4 = which(rcdf >= 0.75, 1)
    rci[1, i] = rgrid[id1[1]]; rci[2, i] = rgrid[id2[1]]
    rci[3, i] = rgrid[id3[1]]; rci[4, i] = rgrid[id4[1]]
  }
  
  # Main outputs: estimates of r and states
  out = list(rmed, rci, rmean, pr, prup, pstate)
  names(out) = c('rmed', 'rci', 'rmean', 'pr', 'prup', 'pstate')
  epiFilterGrowth = out
}