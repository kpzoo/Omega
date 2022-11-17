######################################################################
## Plot EpiFilter estimates and predictions
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” Biorxiv.
######################################################################

# Assumptions
# - uses posterior outputs from either EpiFilter or EpiSmoother

# Inputs - growth rate estimate (rhat) and confidence intervals (rhatci), 
# incidence estimate (Ihat) and confidence intervals (Ihatci), string for figure (plotname),
# true incidenct to compare to (Iplt) and folder path to store results (folres), noise term (eta)

# Output - .eps plots of best r estimate and one-step-ahead incidence predictions

plotEpiFilterGrowth <- function(rhat, rhatci, Inexhat, Inexhatci, plotname, Iplt, folres, eta){
  
  # Check lengths
  if (length(rhat) != length(Inexhat)){
    print(c('rhat length', length(rhat)))
    print(c('Ihat length', length(Inexhat)))
    stop('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    tset = 1:length(rhat)
    
    # Two panel plot of estimates and predictions
    pdf(file=paste0(folres, plotname, '.pdf')) 
    par(mfrow=c(2,1))
    # reprod. num estimates and confidence interval
    plot(tset, rhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (eta = ", eta, ")"), ylab = "reprod. number")
    polygon(c(tset, rev(tset)), c(rhatci[1,], rev(rhatci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(rhatci[3,], rev(rhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    
    # Incidence predictions and confidence interval
    plot(tset, Inexhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (eta = ", eta, ")"), ylab = "incidence", ylim = c(0, max(Iplt)+30))
    polygon(c(tset, rev(tset)), c(Inexhatci[1,], rev(Inexhatci[2,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Inexhatci[3,], rev(Inexhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    points(tset, Iplt, pch = 19, col = 'gray')
    dev.off()
    
  }
}