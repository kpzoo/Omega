#################################################################################
## Compute effective (R) and angular reproduction numbers (Omega) on data 
# from Ali et al 2020 involving changing serial intervals (replicate Fig 2A-2D)
#################################################################################

# Functionality shown in this vignette
# - load COVID daily incidence curves by onset from Ali et al
# - load serial intervals (SIs) over 3 periods from Ali et al
# - estimate R and Omega using EpiFilter (fixed and varying SIs)
# - assumes data is not delayed and any under-reporting is constant

# Clean the workspace and console
closeAllConnections(); rm(list=ls()); cat("\014")
graphics.off(); start_time = as.numeric(Sys.time())

# For fitting serial intervals
library(fitdistrplus)

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/ali/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

# Load incidence data by onset date
casedata = read.csv('data/ali/cases_by_infectee_onset_date.csv')
# Incidence and dates
Iday = casedata$n; dates = casedata$infectee_onsetDate
# Time series lengths
nday = length(dates); tday = 1:nday
# Dates with correction of an error
tdates = as.Date(dates); tdates[length(tdates)] = "2020-02-18"

# Setup grid [Rmin Rmax] noise eta 
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m)

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
# Maximum of prediction grid and confidence (1-conf)
Imax = 2*max(Iday); conf = 0.025

######################################################################
## Smoothed R number estimates under fixed serial interval choices
######################################################################

# Possible fixed serial intervals from Ali et al 2020
nCh = 3; meanCh = c(7.1, 5.2, 3.0); sdCh = c(5.3, 4.7, 4.1)
# Store outputs
Rmeanfix = list(); Rcifix = list(); pRfix1 = list()

# Compute Lday for each serial interval and estimate R
for (j in 1:3) {
  # Specify the fixed serial interval
  pms = c(0,0); mean_si = meanCh[j]; sd_si = sdCh[j]
  pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si
  
  # Approximate serial interval distribution from Ferguson et al
  tdist = c(0, tday); wdist = rep(0, nday)
  for (i in 1:nday){
    wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) - 
      pgamma(tdist[i], shape = pms[1], scale = pms[2])
  }
  
  # Total infectiousness for given serial interval
  Lday = rep(0, nday) 
  for(i in 2:nday){
    # Total infectiousness
    Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
  }
  
  # Filtered (causal) estimates as list [Rmed, Rci, Rmean, pR, pRup, pstate]
  Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], conf)
  # Causal predictions from filtered estimates [pred predci]
  Ifilt = recursPredict(Rgrid, Rfilt$pR, Lday[tday], Rfilt$Rmean, conf, Iday, Imax)
  
  # Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
  Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, nday, Rfilt$pstate, conf)
  # Smoothed predictions from filtered estimates [pred predci]
  Ismooth = recursPredict(Rgrid, Rsmooth$qR, Lday[tday], Rsmooth$Rmean, conf, Iday, Imax)
  
  # Save estimates for plotting
  Rmeanfix[[j]] = Rsmooth$Rmean; Rcifix[[j]] = Rsmooth$Rci 
  
  # Prob of R > 1 (resurgence)
  pRfix = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
  for (i in 1:nday){
    pRfix[i] = sum(Rsmooth$qR[i, id1:m])
  }
  pRfix1[[j]] = pRfix
  # Ouput csv with estimates as well as store to plot
  nam1 = paste0(c(folres, "Rmeanfix", j), collapse = '')
  nam2 = paste0(c(folres, "Rcifix", j), collapse = '')
  write.table(Rsmooth$Rmean, file = nam1, row.names=FALSE, col.names=FALSE)
  write.table(Rsmooth$Rci, file = nam2, row.names=FALSE, col.names=FALSE)
}

######################################################################
## Compute smoothed R estimates with varying serial intervals
######################################################################

# Serial interval changes in 3 periods (pre, post and peak)
sidata = read.csv('data/ali/serial_intervals.csv')
# Unique dates except last one give periods of interest
tuniq = unique(sidata$date_set); tuniq = tuniq[1:3]

# Indices with serial data for each period
si_period = list(); si_id = list(); fit = list()
for (j in 1:3) {
  # Extract relevant serial intervals
  si_id[[j]] = which(sidata$date_set == tuniq[j])
  si_period[[j]] = sidata$si[si_id[[j]]]
  # Fit distributions (matches paper)
  fitj = fitdistr(si_period[[j]], "normal"); fit[[j]] = fitj$estimate
}

# Manual configuring of end points (last extended from 13)
endPt = c("2020/01/22", "2020/01/30", "2020/02/18")

# Apply gamma distributions (even though not true)
wvary = list(); idvary = rep(0, 3)
for (j in 1:3) {
  # Specify the fixed serial interval
  pms = c(0,0); mean_si = fit[[j]][1]; sd_si = fit[[j]][2]
  pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si
  
  # Approximate serial interval distribution from Ferguson et al
  tdist = c(0, tday); wdist = rep(0, nday)
  for (i in 1:nday){
    wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) - 
      pgamma(tdist[i], shape = pms[1], scale = pms[2])
  }
  # Assign this distribution to a period and get end point
  wvary[[j]] = wdist; idvary[j] = which(tdates == endPt[j])
}

# Obtain time varying total infectiousness 
Lvary = rep(0, nday); iset = 1
for(i in 2:nday){
  # Total infectiousness changing with endpoints
  idx = which(i <= idvary); idx = idx[1]
  Lvary[i] = sum(Iday[seq(i-1, 1, -1)]*wvary[[idx]][1:(i-1)])    
}

# Filtered (causal) estimates from time varying serial interval
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lvary[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt$pR, Lvary[tday], Rfilt$Rmean, conf, Iday, Imax)

# Smoothed estimates from time varying serial interval
Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, nday, Rfilt$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth$qR, Lvary[tday], Rsmooth$Rmean, conf, Iday, Imax)

# Prob of R > 1 (resurgence)
pR1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pR1[i] = sum(Rsmooth$qR[i, id1:m])
}

# Save Omega estimates
nam1 = paste0(c(folres, "Rmean"), collapse = '')
nam2 = paste0(c(folres, "Rci"), collapse = '')
write.table(Rsmooth$Rmean, file = nam1, row.names=FALSE, col.names=FALSE)
write.table(Rsmooth$Rci, file = nam2, row.names=FALSE, col.names=FALSE)

######################################################################
## Compute smoothed Omega estimates which do not need serial intervals
######################################################################

# Root mean square incidence for Omega under window delta
Irms = rep(0, nday); delta = 15
for(i in 2:nday){
  # Edge effect cases
  if(i-1 < delta){
    Irms[i] = sqrt(sum(Iday[seq(i-1, 1, -1)]^2))/sqrt(i-1)     
  }else{
    # Omega only truly valid from this point
    Irms[i] = sqrt(sum(Iday[seq(i-1, i-delta, -1)]^2))/sqrt(delta) 
  }
}

# Filtered (causal) estimates as list [Ommed, Omci, Ommean, pOm, pOmup, pOmstate]
Omfilt = epiFilter(Rgrid, m, eta, pR0, nday, Irms[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
IfiltOm = recursPredict(Rgrid, Omfilt$pR, Irms[tday], Omfilt$Rmean, conf, Iday, Imax)

# Smoothed estimates as list of [Ommed, Omhatci, Ommean, qOm]
Omsmooth = epiSmoother(Rgrid, m, Omfilt$pR, Omfilt$pRup, nday, Omfilt$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
IsmoothOm = recursPredict(Rgrid, Omsmooth$qR, Irms[tday], Omsmooth$Rmean, conf, Iday, Imax)

# Prob of R > 1 (resurgence)
pOm1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pOm1[i] = sum(Omsmooth$qR[i, id1:m])
}

# Save Omega estimates
nam1 = paste0(c(folres, "Ommean"), collapse = '')
nam2 = paste0(c(folres, "Omci"), collapse = '')
write.table(Omsmooth$Rmean, file = nam1, row.names=FALSE, col.names=FALSE)
write.table(Omsmooth$Rci, file = nam2, row.names=FALSE, col.names=FALSE)

######################################################################
## Visualisation and ouput of results (as csv files)
######################################################################

# Plotting to match time frame of Ali et al
idplt = 5:(nday-5); tset = tdates[idplt]
# Adjust SI change time indices
idchplt = rep(0, 2)
for (j in 1:2) {
  idchplt[j] = which(tset == endPt[j])
}

Rmeanfixplt = list(); Rcifixplt = list(); pRfix1plt = list()
for (j in 1:3) {
  Rmeanfixplt[[j]] = Rmeanfix[[j]][idplt]; Rcifixplt[[j]] = Rcifix[[j]][, idplt]
  pRfix1plt[[j]] = pRfix1[[j]][idplt]
}

plotOmegaCaseStudy(Rsmooth$Rmean[idplt], Rsmooth$Rci[, idplt], Rmeanfixplt, Rcifixplt, idchplt,
                   Omsmooth$Rmean[idplt], Omsmooth$Rci[, idplt], pRfix1plt, pR1[idplt], pOm1[idplt],
                   meanCh, tset, 'compareAll', folres)

# Log running time in minutes
end_time = as.numeric(Sys.time()); exec_time = end_time - start_time 
exec_time = exec_time/60; exec_time = round(exec_time, 3)
print(paste0(c("Completed in ", exec_time, " mins"), collapse = ""))