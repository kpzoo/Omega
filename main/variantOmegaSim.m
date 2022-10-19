% Simulate two epidemics and estimate metrics
function [simOut, estOut] = variantOmegaSim(epiNo, scenNo, nday, I0, idcutChoice, simVals, delw)

% Assumptions and notes
% - examine true values and estimates (EpiFilter)
% - change R and generation time distribution
% - decide if including predictions or not
% - only setup for nw = 2 variants

% Two generation time distributions with different means
nw = 2; meanPmFrac = [1 delw]; tday = 1:nday;
wch = zeros(nw, nday); idcut = zeros(1, nw); Rtrue = wch;

% Obtain generation time distributions and true R
for i = 1:nw
    wch(i, :) = getGenTimeDistVary(nday, epiNo, meanPmFrac(i));
    % Find cutoff of generation time distribution
    idcut(i) = find(abs(cumsum(wch(i, :))-1) < 10^-8, 1, 'first');
    % Define a scenario to simulate (true R)
    Rtrue(i, :) = scenarioReprod(scenNo, nday, simVals{i});
end

% Means and variances of distributions
wmeans = wch*tday'; wmeans = wmeans';
wvars = (wch*(tday'.^2))' - wmeans.^2;

% Generation time properties disease based (only EVD and COVID)
switch(epiNo)
    case 4
        % COVID-19 growth parameters
        shapePm = (wmeans.^2)/(4.23^2); scalePm = wmeans./shapePm;
    case 5
        % EVD growth parameters
        shapePm = (wmeans.^2)/(9.3^2); scalePm = wmeans./shapePm;
end

% Case counts, infectiousness and initial value
Iday = zeros(nw, nday); Iday(:, 1) = I0; 
% Norms, incidence 2-mean, omega, growth rate, infectiousness
wnorm = zeros(nw, nday); rtrue = wnorm; omega = wnorm; 
Imax = zeros(1, nw); Ipmean = Iday; Lday = Iday; Imean = Imax;


% Simulate epidemics and compute statistics
for j = 1:nw
    while Imax(j) < 200 || Imean(j) < 50
        for i = 2:nday
            % Relevant part of variant generation times
            if i-1 <= idcut(j)
                wt = wch(j, 1:i-1); It = Iday(j, i-1:-1:1);
            else
                wt = wch(j, 1:idcut(j));
                It = Iday(j, i-1:-1:(i-1-idcut(j)+1));
            end

            % Relevant part of generation time for omega estimation
            if i-1 <= idcutChoice
                Itnorm = Iday(j, i-1:-1:1);
            else
                Itnorm = Iday(j, i-1:-1:(i-1-idcutChoice+1));
            end

            % Total infectiousness and 2-mean denominators
            Lday(j, i) = sum(It.*wt);
            Ipmean(j, i) = norm(Itnorm)/sqrt(min(i-1, idcutChoice));

            % Renewal incidence simulation draw
            Iday(j, i) = poissrnd(Rtrue(j, i)*Lday(j, i));
            
            % True omega metric and growth rate
            omega(j, i) = (Rtrue(j, i)*Lday(j, i))/Ipmean(j, i);
            rtrue(j, i) = (Rtrue(j, i).^(1/shapePm(j)) - 1)/scalePm(j);
        end
        % Maximum of simulation
         Imax(j) = max(Iday(j, :)); disp(['Imax = ' num2str(Imax(j))]);
         % Mean of simulation
         Imean(j) = mean(Iday(j, :)); disp(['Imean = ' num2str(Imean(j))]);
    end
end

% Output key simulation variables
simOut.Iday = Iday; simOut.Lday = Lday;
simOut.Ipmean = Ipmean; simOut.tday = tday;
simOut.Rtrue = Rtrue; simOut.rtrue = rtrue; 
simOut.wch = wch; simOut.wmeans = wmeans;
simOut.wvars = wvars; simOut.omega = omega;
simOut.meanPmFrac = meanPmFrac; simOut.Imax = Imax;

% Estimate metrics for variants 

% Grid limits and noise level for R estimation
Rmin = 0.01; Rmax = 10; eta = 0.1; 
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Estimates from each epidemic
Rlow = zeros(nw, nday); Rhigh = Rlow; Rmean = Rlow; 
qR = cell(1, nw); Omlow = Rlow; Omhigh = Rlow; Ommean = Rlow; 
qOm = qR; rmean = Rlow; rlow = Rlow; rhigh = Rlow; 

% Estimates of Omega, growth rates and R
for j = 1:nw
    % EpiFilter smoothing and posterior prediction
    [~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lday(j, :), Iday(j, :));
    [~, Rlow(j, :), Rhigh(j, :), Rmean(j, :), qR{j}] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

    % EpiFilter omega estimation
    [~, ~, ~, ~, pOm, pOmup, pOmstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Ipmean(j, :), Iday(j, :));
    [~, Omlow(j, :), Omhigh(j, :), Ommean(j, :), qOm{j}] = runEpiSmoother(Rgrid, m, nday, pOm, pOmup, pOmstate);

    % Growth rate estimates
    rmean(j, :) = (Rmean(j, :).^(1/shapePm(j)) - 1)/scalePm(j);
    rlow(j, :) = (Rlow(j, :).^(1/shapePm(j)) - 1)/scalePm(j);
    rhigh(j, :) = (Rhigh(j, :).^(1/shapePm(j)) - 1)/scalePm(j);
end

% Output key estimation variables
estOut.Rmean = Rmean; estOut.Rlow = Rlow; estOut.Rhigh = Rhigh;
estOut.rmean = rmean; estOut.rlow = rlow; estOut.rhigh = rhigh;
estOut.Ommean = Ommean; estOut.Omlow = Omlow; estOut.Omhigh = Omhigh;




