######################################################################
## Plot EpiFilter estimates and predictions for R and Omega
######################################################################

# Assumptions
# - uses posterior outputs from either EpiFilter or EpiSmoother

# Inputs - reproduction number estimate (Rhat) and confidence intervals (Rhatci), 
# incidence estimate (Ihat) and confidence intervals (Ihatci), string for figure (plotname),
# true incidenct to compare to (Iplt) and folder path to store results (folres), noise term (eta)

# Output - .eps plots of best R estimate and one-step-ahead incidence predictions

plotOmegaR <- function(Rhat, Rhatci, IRhat, IRhatci, Omhat, Omhatci, IOmhat, IOmhatci,
                          plotname, Iplt, folres, eta, delta){
  
  # Check lengths
  if (length(Rhat) != length(IRhat) | length(Omhat) != length(IOmhat)){
    print(c('[Rhat Omhat] lengths', length(Rhat), length(Omhat)))
    print(c('[Ihat Iomhat] lengths', length(IRhat), length(IOmhat)))
    stop('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    tset = 1:length(Rhat)
    
    # Two panel plot of estimates and predictions
    pdf(file=paste0(folres, plotname, '.pdf')) 
    par(mfrow=c(2,1))
    
    # Reprod. num estimates and confidence interval
    plot(tset, Rhat, type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, max(Rhatci[2,]) + 0.2),
         xlab = paste0("time (eta = ", eta, ")"), ylab = "reprod. number")
    polygon(c(tset, rev(tset)), c(Rhatci[1,], rev(Rhatci[2,])), 
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rhatci[3,], rev(Rhatci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    
    # Omega estimates and confidence interval
    lines(tset, Omhat, type = 'l', bty = 'l', lwd = 2, col='blue')
    polygon(c(tset, rev(tset)), c(Omhatci[1,], rev(Omhatci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Omhatci[3,], rev(Omhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    lines(c(delta, delta), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    
    # Incidence predictions and confidence interval
    plot(tset, IRhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (eta = ", eta, ")"), ylab = "incidence", ylim = c(0, max(Iplt)+30))
    polygon(c(tset, rev(tset)), c(IRhatci[1,], rev(IRhatci[2,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(IRhatci[3,], rev(IRhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    
    lines(tset, IRhat, type = 'l', bty = 'l', lwd = 2, col='red')
    polygon(c(tset, rev(tset)), c(IRhatci[1,], rev(IRhatci[2,])),
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(IRhatci[3,], rev(IRhatci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    
    points(tset, Iplt, pch = 19, col = 'gray')
    dev.off()
    
    # Single plot of estimates 
    pdf(file=paste0(folres, paste0(c(plotname, "_est"), collapse = ""), '.pdf')) 
    par(mfrow=c(1,1))
    
    # Reprod. num estimates and confidence interval
    plot(tset, Rhat, type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, max(Rhatci[2,]) + 0.2),
         xlab = "time (days)", ylab = "reprod. number, omega")
    polygon(c(tset, rev(tset)), c(Rhatci[1,], rev(Rhatci[2,])), 
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rhatci[3,], rev(Rhatci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    
    # Omega estimates and confidence interval
    lines(tset, Omhat, type = 'l', bty = 'l', lwd = 2, col='blue')
    polygon(c(tset, rev(tset)), c(Omhatci[1,], rev(Omhatci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Omhatci[3,], rev(Omhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    lines(c(delta, delta), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    dev.off()
    
  }
}