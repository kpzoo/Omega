#################################################################################
## Compute percent confidence at any time that an epidemic is eliminated or over
# From: Parag KV, Cowling BJ, Donnelly CA (2021). Deciphering early-warning signals 
# of SARS-CoV-2 elimination and resurgence from limited data at multiple scales. 
# J. R. Soc. Interface 18(185):20210444..
#################################################################################

# Notes on usage
# - assumes some gamma shaped serial interval but can be generalised
# - uses diffusion-Poisson assumptions of EpiFilter
# - if no local-import distinctions then set Iloc = Iday

# Inputs - total incidence (Iday), local incidence (Iloc), grid on R space (Rgrid),
# dimension of grid (m), shape-scale parameters of gamma serial interval (pms)
# state noise (eta), prior on R (pR0) and confidence (1-conf) for R estimates

# Output - confidence in elimination at end-time of input incidence curve under
# the filtered R (z0) and the smoothed R (z1) in %

computeZ <- function(Iloc, Iday, Rgrid, m, pms, eta, pR0, conf){
  
  # Test input dimensions
  if(length(Iloc) != length(Iday) | length(pR0) != m | length(Rgrid) != m){
    stop("Incorrect input dimensions")
  }
  
  # Define a zero look-ahead sequence
  Iz = rep(0, 200); lz = length(Iz) - 1
  
  # Append epidemic curve with pseudo-data
  Icurr = c(Iloc, Iz); ncurr = length(Icurr)
  # Total cases (local and imported) for computing Lam
  Ilamcurr = c(Iday, Iz); tcurr = 1:ncurr
  
  # Range of index time points
  idz = length(Iloc) + 1; ir = idz:(ncurr-1) 
  
  # Gamma serial interval (should match that used in EpiFilter)
  tdist = c(0, tcurr); wdist = rep(0, nday)
  for (i in 1:ncurr){
    wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) - 
      pgamma(tdist[i], shape = pms[1], scale = pms[2])
  }
  
  # Compute successive total infectiousness from total cases
  Lcurr = rep(0, ncurr);
  for(i in 2:ncurr){
    # Relevant part of SI: Pomega(1:i-1))
    Lcurr[i] = sum(Ilamcurr[seq(i-1, 1, by = -1)]*wdist[1:(i-1)]);
  }
  
  # Run EpiFilter and smoother on relevant incidence and total infectiousness
  Rfilt = epiFilter(Rgrid, m, eta, pR0, ncurr, Lcurr, Icurr, conf)
  Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, ncurr, Rfilt$pstate, conf)
  # Filtered and smoothed posterior distributions (pR and qR)
  pR = Rfilt$pR; qR = Rsmooth$qR
  
  # Sequences of probabilities across zeros
  zseq0 = rep(0, lz); zseq1 = zseq0; ix = 1;
  for(i in ir){
    # Lam and elimin prob at this time
    L = Lcurr[i+1]; zL = exp(-L*Rgrid)
    
    # Integrate over grid for filtered Z number sequence
    zseq0[ix] = sum(zL*pR[i, ])
    # Integrate over grid for smoothed Z number sequence
    zseq1[ix] = sum(zL*qR[i, ])
    ix = ix + 1
  }

  # Estimated probability of elimination at this time
  z0 = prod(zseq0); z1 = prod(zseq1);
  # Main outputs: filtered and smoothed confidence scores
  computeZ = c(100*z0, 100*z1)
}

