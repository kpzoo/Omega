% Simulate single epidemic and estimate metrics
function [simOut, estOut] = singleOmegaSim(epiNo, scenNo, nday, I0, idcutChoice, gendel, deltype)

% Assumptions and notes
% - examine true values and estimates (EpiFilter)
% - change R and generation time distribution
% - decide if including predictions or not

% Case counts and initial value
Iday = zeros(1, nday); Iday(1) = I0; 
% Total infectiousness and time of runs
Lday = Iday; tday = 1:nday;

% Define a scenario to simulate (true R)
simVals = setupScenario(scenNo);
Rtrue = scenarioReprod(scenNo, nday, simVals);

% Define a number of generation time fractions
if deltype == 0
    % Alternating changes
    meanPmFrac = [1 1-gendel 1-gendel 1+gendel 1-gendel];
    nwChoice = 4; nw = length(meanPmFrac);
else
    % Persistent changes
    meanPmFrac = [1 1 1 1+gendel 1+gendel];
    meanPmFrac = [1 1 1 1-gendel 1-gendel];
    nwChoice = 1; nw = length(meanPmFrac);
end

% Uniformly spread times for initiating these fractions
wchtime = round(linspace(1, nday, nw+1)); wchtime = wchtime(1:end-1);

% Obtain these generation time distributions
wch = zeros(nw, nday); idcut = zeros(1, nw);
for i = 1:nw
    wch(i, :) = getGenTimeDistVary(nday, epiNo, meanPmFrac(i));
    % Find cutoff of generation time distribution 
    idcut(i) = find(abs(cumsum(wch(i, :))-1) < 10^-8, 1, 'first');
end

% Means of distributions
wmeans = wch*tday'; wmeans = wmeans';
% Variances of distributions
wvars = (wch*(tday'.^2))' - wmeans.^2;

% Norms on generation time and incidence
wnorm = zeros(1, nday); omega = wnorm; Ipmean = wnorm; 

% Initialise changes and keep original distribution and Lam
wtorig = wch(1, :); Lorig = Lday; Imax = 0;

% Simulate epidemics and compute statistics
while Imax < 200
    id = 1;
    for i = 2:nday
        % Choose generation time distribution
        if any(i == wchtime)
            id = id + 1;
        end

        % Relevant part of altered generation time distributions
        if i-1 <= idcut(id)
            wt = wch(id, 1:i-1); It = Iday(i-1:-1:1);
        else
            wt = wch(id, 1:idcut(id));
            It = Iday(i-1:-1:(i-1-idcut(id)+1));
        end

        % Relevant part of generation time for omega estimation
        if i-1 <= idcutChoice
            Itnorm = Iday(i-1:-1:1);
        else
            Itnorm = Iday(i-1:-1:(i-1-idcutChoice+1));
        end
        Ipmean(i) = norm(Itnorm)/sqrt(min(i-1, idcutChoice));

        % Total infectiousness (true and original)
        Lday(i) = sum(It.*wt);
        Lorig(i) = sum(Iday(i-1:-1:1).*wtorig(1:i-1));

        % Renewal incidence and power mean
        Iday(i) = poissrnd(Rtrue(i)*Lday(i));

        % Left and right side of omega metric
        omega(i) = (Rtrue(i)*Lday(i))/Ipmean(i);
    end
    % Maximum of simulation
    Imax = max(Iday); disp(['Imax = ' num2str(Imax)]);
end

% Output key simulation variables
simOut.Iday = Iday; simOut.Lday = Lday;
simOut.Lorig = Lorig; simOut.Ipmean = Ipmean;
simOut.Rtrue = Rtrue; simOut.tday = tday;
simOut.wch = wch; simOut.wmeans = wmeans;
simOut.wvars = wvars; simOut.wchtime = wchtime;
simOut.meanPmFrac = meanPmFrac; simOut.omega = omega;
simOut.nwChoice = nwChoice; simOut.Imax = Imax;

% Estimate metrics for fixed generation times

% Grid limits and noise level for R estimation
Rmin = 0.01; Rmax = 10; eta = 0.1; 

% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter smoothing and posterior prediction
[~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lorig, Iday);
[~, Rlow, Rhigh, Rmean, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiFilter omega estimation
[~, ~, ~, ~, pOm, pOmup, pOmstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Ipmean, Iday);
[~, Omlow, Omhigh, Ommean, qOm] = runEpiSmoother(Rgrid, m, nday, pOm, pOmup, pOmstate);


% Complementary CDFs for R > 1
id1 = find(Rgrid <= 1, 1, 'last');
Fq(:, i) = 1 - sum(qR(:, 1:id1), 2)'; FqOm(:, i) = 1 - sum(qOm(:, 1:id1), 2)'; 

% Mean and quantiles
Fqm = mean(Fq, 2); FqmOm = mean(FqOm, 2); 
Fqq = quantile(Fq, [0.025 0.975], 2)'; FqqOm = quantile(FqOm, [0.025 0.975], 2)'; 

% Output key estimation variables
estOut.Rmean = Rmean; estOut.Ommean = Ommean;
estOut.Rlow = Rlow; estOut.Rhigh = Rhigh;
estOut.Omlow = Omlow; estOut.Omhigh = Omhigh;
estOut.Fq = Fq; estOut.Fqm = Fqm; estOut.Fqq = Fqq;
estOut.FqOm = FqOm; estOut.FqmOm = FqmOm; estOut.FqqOm = FqqOm;

