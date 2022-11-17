######################################################################
## Plot EpiFilter estimates for R and Omega on Ali et al data
######################################################################

# Assumptions
# - uses posterior outputs from EpiFilter and assumes 3 fixed and 1 varying SI based R
# Inputs - reproduction number estimate (Rmean) and confidence intervals (Rci), 
# estimate (Rmeanfix) and confidence intervals (Rcifix) under fixed SIs, 
# estimate (Ommean) and confidence intervals (Omci) of Omega metric, 
# string for figure (plotname) and folder path to store results (folres)

# Output - .eps plots of best R and Omega estimates

plotOmegaCaseStudy <- function(Rmean, Rci, Rmeanfix, Rcifix, idchplt, Ommean, Omci, pRfix1, pR1, pOm1,
                               meanSI, tset, plotname, folres){
  
  # Check lengths
  if (length(Rmean) != length(Rmeanfix[[1]]) | length(Ommean) != length(Rmean)){
    print(c('[Rmean Ommean] lengths', length(Rmean), length(Ommean)))
    print(c('[Rmean Rmeanfix] lengths', length(Rmean), length(Rmeanfix[[1]])))
    stop('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    #tset = 1:length(Rmean)
    
    # Four panel plot of estimates against Omega
    pdf(file=paste0(folres, plotname, '.pdf')) 
    par(mfrow=c(2,2))
    
    # Reprod. num estimates under assumed fixed serial intervals
    for (j in 1:3) {
      plot(tset, Rmeanfix[[j]], type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, max(Rcifix[[j]][2,]) + 0.2),
           xlab = paste0("time (g = ", meanSI[j], ")"), ylab = "reprod. number")
      polygon(c(tset, rev(tset)), c(Rcifix[[j]][1,], rev(Rcifix[[j]][2,])), 
              col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
      polygon(c(tset, rev(tset)), c(Rcifix[[j]][3,], rev(Rcifix[[j]][4,])), 
              col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
      
      # Omega estimates and confidence interval
      lines(tset, Ommean, type = 'l', bty = 'l', lwd = 2, col='blue')
      polygon(c(tset, rev(tset)), c(Omci[1,], rev(Omci[2,])), 
              col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
      polygon(c(tset, rev(tset)), c(Omci[3,], rev(Omci[4,])), 
              col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
      lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
      lines(c(tset[idchplt[1]], tset[idchplt[1]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
      lines(c(tset[idchplt[2]], tset[idchplt[2]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    }
    
    # Reprod. num estimates under time-varying serial intervals
    plot(tset, Rmean, type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, max(Rci[2,]) + 0.2),
         xlab = paste0("time (g = ", "variable", ")"), ylab = "reprod. number")
    polygon(c(tset, rev(tset)), c(Rci[1,], rev(Rci[2,])), 
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rci[3,], rev(Rci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    
    # Omega estimates and confidence interval
    lines(tset, Ommean, type = 'l', bty = 'l', lwd = 2, col='blue')
    polygon(c(tset, rev(tset)), c(Omci[1,], rev(Omci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Omci[3,], rev(Omci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    lines(c(tset[idchplt[1]], tset[idchplt[1]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    lines(c(tset[idchplt[2]], tset[idchplt[2]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    dev.off()
    
    # Single plot of estimates 
    pdf(file=paste0(folres, paste0(c(plotname, "_est"), collapse = ""), '.pdf')) 
    par(mfrow=c(1,1))
    
    # Reprod. num estimates and confidence interval
    plot(tset, Rmean, type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, max(Rci[2,]) + 0.2),
         xlab = "time (days)", ylab = "reprod. number, omega")
    polygon(c(tset, rev(tset)), c(Rci[1,], rev(Rci[2,])), 
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rci[3,], rev(Rci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    
    # Omega estimates and confidence interval
    lines(tset, Ommean, type = 'l', bty = 'l', lwd = 2, col='blue')
    polygon(c(tset, rev(tset)), c(Omci[1,], rev(Omci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Omci[3,], rev(Omci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    lines(c(tset[idchplt[1]], tset[idchplt[1]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    lines(c(tset[idchplt[2]], tset[idchplt[2]]), c(0, 3), lwd = 2, col = 'lightgrey', lty = 'dashed')
    dev.off()
    
    # Single plot of p(R > 1)
    pdf(file=paste0(folres, paste0(c(plotname, "_R1"), collapse = ""), '.pdf')) 
    par(mfrow=c(1,1))
    plot(tset, pR1, type = 'l', bty = 'l', lwd = 2, col='red', ylim = c(0, 1.05),
         xlab = "time (days)", ylab = "reprod. number, omega")
    for (j in 1:3) {
      lines(tset, pRfix1[[j]], lwd = 2, col = 'lightgrey')
    }
    lines(tset, pOm1, lwd = 2, col = 'blue')
    dev.off()
  }
}