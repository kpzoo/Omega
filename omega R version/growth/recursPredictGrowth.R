######################################################################
## Bayesian recursive prediction via EpiFilter
# From: Parag, KV (2021). Improved real-time estimation of reproduction numbers at
# low case incidence and between epidemic waves. PLOS Comput Biol 17(9):e1009347.
######################################################################

# Notes and assumptions
# - observation model is Poisson over the growth rate and past incidence
# - growth rate state space model is a simple diffusion
# - can apply causal or smoothing posteriors over r to predict incidence
# - sets maximum on possible Igrid to interrogate as Imax

# Inputs - grid on reproduction numbers (rgrid), posterior on r (pr), 
# incidence (Iday), mean r esimate using pr (rmean), confidence level (a and 50%)

# Output - mean prediction (pred), confidence intervals (predci), predictive score (score)

recursPredictGrowth <- function(rgrid, pr, rmean, a, Iday, Imax){
  
  # Grid size and length of time series
  nday = nrow(pr); m = ncol(pr)
  # Test lengths of inputs
  if (length(rgrid) != m | length(Iday) != nday){
    stop("Input vectors of incorrect dimension")
  }
  
  # Mean prediction: Iday[i] => Iday[i+1]
  pred = Iday*exp(rmean); pred = pred[1:length(pred)-1]
  
  # Discrete space of possible predictions
  Igrid = 0:Imax; lenI = length(Igrid);
  
  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop("Epidemic size too large")  
  }
  
  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predci = matrix(0, 4, nday-1)
  # Score of fit and predictive distributions
  score = rep(0, nday-1); pIdist = matrix(0, nday-1, lenI)
  
  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Iday[i]*exp(rgrid)
    # Prob of any I marginalised over rgrid
    pI = rep(0, lenI)
    
    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of r
      pI[j] = sum(pIset*pr[i, ])
    }
    pIdist[i,] = pI
    
    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)
    
    # Score of predictions
    score[i] = (Iday[i+1] - pred[i])^2
    #idI = which(Igrid == Iday[i+1]); pid = pI[idI]
    #if(pid == 0){score[i] = 0}else{score[i] = -log(pid)}
    
    # Assign prediction results
    predci[1, i] = Igrid[id1[1]]; predci[2, i] = Igrid[id2[1]]
    predci[3, i] = Igrid[id3[1]]; predci[4, i] = Igrid[id4[1]]
  }
  
  # Main outputs: mean and 50% and (1-a)% predictions
  out = list(pred, predci, score, pIdist)
  names(out) = c('pred', 'predci', 'score', 'dist')
  recursPredictGrowth = out
}