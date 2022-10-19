% Estimate metrics from Ali et al data and serial intervals
function [simOut, estOut] = hkOmega(nday, Iday, delta, wmeans, wvars, idch)

% Assumptions and notes
% - compute R and omega (window delta) using EpiFilter
% - changing serial interval distribution (idch is index)

% Total infectiousness and time of runs
Lday = Iday; tday = 1:nday;

% Generation time distribution statistics
nw = length(idch); wch = zeros(nw, nday); 
for i = 1:nw
    % Shape and scale of gamma distribution
    scalePm = wvars(i)/wmeans(i); shapePm = wmeans(i)/scalePm;
    wch(i, :) = gamcdf(tday, shapePm, scalePm) - gamcdf(tday-1, shapePm, scalePm);
end

% Compute denominators for changing generation times
id = 1; Ipmean = Lday; Lfix = zeros(nw, nday); Lfix(:, 1) = Lday(1);
for i = 2:nday
    % Choose generation time distribution
    if any(i == idch+1)
        id = id + 1;
    end

    % RMS of infections for omega given window delta
    if i-1 <= delta
        Itnorm = Iday(i-1:-1:1);
    else
        Itnorm = Iday(i-1:-1:(i-1-delta+1));
    end
    Ipmean(i) = norm(Itnorm)/sqrt(min(i-1, delta));

    % Total infectiousness from data for R
    Lday(i) = sum(Iday(i-1:-1:1).*wch(id, 1:i-1));
end

% Total infectiousness assuming one of the fixed distributions
for j = 1:nw
    for i = 2:nday
        Lfix(j, i) = sum(Iday(i-1:-1:1).*wch(j, 1:i-1));
    end
end

% Output key simulation variables
simOut.Iday = Iday; simOut.Lday = Lday; simOut.Lfix = Lfix;
simOut.Ipmean = Ipmean; simOut.tday = tday; simOut.wch = wch;
simOut.wmeans = wmeans; simOut.wvars = wvars; simOut.wchtime = idch;

%% Estimate metrics for fixed and varying generation times

% Grid limits and noise level for R estimation
Rmin = 0.01; Rmax = 10; eta = 0.1; 
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter for fixed generation times
Rfix = zeros(nw, nday); Rlfix = Rfix; Rhfix = Rfix;
for i = 1:nw
    [~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lfix(i, :), Iday);
    [~, Rlfix(i, :), Rhfix(i, :), Rfix(i, :), ~] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);
end

% EpiFilter for varying generation times
[~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lday, Iday);
[~, Rlow, Rhigh, Rmean, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiFilter omega estimation
[~, ~, ~, ~, pOm, pOmup, pOmstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Ipmean, Iday);
[~, Omlow, Omhigh, Ommean, qOm] = runEpiSmoother(Rgrid, m, nday, pOm, pOmup, pOmstate);

% Output key estimation variables
estOut.Rlow = Rlow; estOut.Rhigh = Rhigh; estOut.Rmean = Rmean; 
estOut.Omlow = Omlow; estOut.Omhigh = Omhigh; estOut.Ommean = Ommean;
estOut.Rlfix = Rlfix; estOut.Rhfix = Rhfix; estOut.Rfix = Rfix;

%% Estimate growth rates and another omega


% Growth rate estimates from R with changing serial interval
rmean = zeros(1, nday); rlow = rmean; rhigh = rmean;

% Compute estimates from R given changetimes
idvar = [0 idch];
for j = 1:nw
    scalePm = wvars(j)/wmeans(j); shapePm = wmeans(j)/scalePm;
    % Time frame for this generation time
    tframe = idvar(j)+1:idvar(j+1);
    % Gamma distributiin based growth rates
    rmean(tframe) = (Rmean(tframe).^(1/shapePm) - 1)/scalePm;
    rlow(tframe) = (Rlow(tframe).^(1/shapePm) - 1)/scalePm;
    rhigh(tframe) = (Rhigh(tframe).^(1/shapePm) - 1)/scalePm;
end

% Compute estimates from omega given delta
A = 1/(2*delta); B = Ommean.^2;
rommean = getGrowthfromOmega(A, B, Ommean);

% Transform the confidence intervals naively
B = Omlow.^2; romlow = getGrowthfromOmega(A, B, Omlow);
B = Omhigh.^2; romhigh = getGrowthfromOmega(A, B, Omhigh);

% Output key estimation variables
estOut.rlow = rlow; estOut.rhigh = rhigh; estOut.rmean = rmean; 
estOut.romlow = romlow; estOut.romhigh = romhigh; estOut.rommean = rommean;



    

