#################################################################################
## Compute effective (R) and angular reproduction numbers (Omega) on simulalted
# data with no variations in generation time
#################################################################################

# Clean the workspace and console
closeAllConnections(); rm(list=ls()); cat("\014")
graphics.off(); start_time = as.numeric(Sys.time())

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/sim/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

library(plotly)

#################################################################################
################ Simulate epidemic data from a sinusoidal R #####################

# Length of simulated time series and true R
nday = 200; tday = 1:nday
# Choose a simulation example
simType = 0
if(simType){
  R = 1.3 + 1.2*sin(4*(pi/180)*tday)
  } else{
  R = rep(1, nday); R[1:60] = 2; R[61:90] = 0.5; R[91:nday] = 1
}

# Shape-scale parameters of gamma serial interval
mean_si <- 6.5; sd_si <- 4.03; pms = c(0,0)
pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si

# Approximate serial interval distribution from Ferguson et al
tdist = c(0, tday); wdist = rep(0, nday)
for (i in 1:nday){
  wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) - 
    pgamma(tdist[i], shape = pms[1], scale = pms[2])
}

# Total infectiousness and incidence
Lday = rep(0, nday); Iday = Lday
# Initialise some infections
Iday[1] = 5; Lday[1] = 5

# Simulate from standard renewal model
for(i in 2:nday){
  # Total infectiousness is a convolution
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])   
  # Poisson renewal model
  Iday[i] = rpois(1, Lday[i]*R[i])
}

# For window delta compute rmse infections
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

#################################################################################
################ Estimate R and Omega using EpiFilter ###########################

# Setup grid [Rmin Rmax] noise eta and CIs conf
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m)

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
# Maximum of prediction grid and confidence (1-conf)
Imax = 2*max(Iday); conf = 0.025

# Filtered (causal) estimates as list [Rmed, Rci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt$pR, Lday[tday], Rfilt$Rmean, conf, Iday, Imax)

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, nday, Rfilt$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth$qR, Lday[tday], Rsmooth$Rmean, conf, Iday, Imax)
# Prob of R > 1 (resurgence)
pR1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pR1[i] = sum(Rsmooth$qR[i, id1:m])
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

#################################################################################
########################### Visualise outputs ###################################


# Plot and save estimates and predictions from smoothing
plotEpiFilterSim(Rsmooth$Rmean[2:nday], Rsmooth$Rci[, 2:nday], Ismooth$pred, Ismooth$predci,
              'RSmooth', Iday[2:nday], folres, eta, R[2:nday])

plotEpiFilterSim(Omsmooth$Rmean[2:nday], Omsmooth$Rci[, 2:nday], IsmoothOm$pred, IsmoothOm$predci,
              'OmSmooth', Iday[2:nday], folres, eta, R[2:nday])

plotOmegaRSim(Rsmooth$Rmean[2:nday], Rsmooth$Rci[, 2:nday], Ismooth$pred, Ismooth$predci,
           Omsmooth$Rmean[2:nday], Omsmooth$Rci[, 2:nday], IsmoothOm$pred, IsmoothOm$predci,
           'compareROmega', Iday[2:nday], folres, eta, delta, R[2:nday], 0)

# Dates for plotting
tdates = tday

# Smoothed local R estimates with confidence intervals
fig = plot_ly(x = tdates, y = Rsmooth$Rci[1,], type = 'scatter', mode = 'lines', name = 'lower CI (R)',
              showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = Rsmooth$Rci[2,], name = 'upper CI (R)', fill = 'tonexty',
                        showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = rep(1, nday), name = 'R = 1', line = list(color = 'lightblue'))
fig = fig %>% add_trace(x = tdates, y = R, name = 'true R', line = list(color = 'black'))
fig = fig %>% layout(xaxis = list(title = 'time (days)' ), yaxis = list(title = 'reproduction number'))
ymax = max(Rsmooth$Rci[2,]) + 0.2; figRest = fig; suppressWarnings(print(figRest))

# Smoothed local Omega estimates with confidence intervals
fig = plot_ly(x = tdates, y = Omsmooth$Rci[1,], type = 'scatter', mode = 'lines', name = 'lower CI (R)',
              showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = Omsmooth$Rci[2,], name = 'upper CI (R)', fill = 'tonexty',
                        showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = rep(1, nday), name = 'R = 1', line = list(color = 'lightblue'))
fig = fig %>% add_trace(x = tdates, y = R, name = 'true R', line = list(color = 'black'))
fig = fig %>% layout(xaxis = list(title = 'time (days)' ), yaxis = list(title = 'angular reproduction number'))
ymax = max(Omsmooth$Rci[2,]) + 0.2; figRest = fig; suppressWarnings(print(figRest))

# Log running time in minutes
end_time = as.numeric(Sys.time()); exec_time = end_time - start_time 
exec_time = exec_time/60; exec_time = round(exec_time, 3)
print(paste0(c("Completed in ", exec_time, " mins"), collapse = ""))
