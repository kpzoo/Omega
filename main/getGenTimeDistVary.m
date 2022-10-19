% Simulate types of generation time distribution (fixed)
function w = getGenTimeDistVary(tmax, epiNo, meanCh)

% Assumptions and notes
% - allows alteration of mean of generation time distributions
% - includes different disease type scenarios
% - gamma for 5 different diseases
% - other forms including uniform and flare (bimodal)

% Possible scenarios available - must match calling function
epiNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'EVD', 'Geometric',...
    'Uniform', 'Bimodal', 'Window'};

%% Generation time distribution specification

% Hyerparameters of serial distribution (must be gamma)
switch(epiNo)
    case 1
        % Marburg distribution from van Kerkhove 2015
        wmean = 9*meanCh; shapePm = (wmean^2)/(5.4^2); 
        w = gammaDistr(1:tmax, wmean, shapePm);
    case 2
        % MERS distribution from Cauchemez 2016
        wmean = 6.8*meanCh; shapePm = (wmean^2)/(4.1^2);
        w = gammaDistr(1:tmax, wmean, shapePm);
    case 3
        % Measles distribution from Cori 2013
        wmean = 14.9*meanCh; shapePm = (wmean^2)/(3.9^2);
        w = gammaDistr(1:tmax, wmean, shapePm);
    case 4
        % COVID-19 distribution from Ferguson 2020
        %wmean = 6.5*meanCh; shapePm = (1/0.65)^2;
        wmean = 6.5*meanCh; shapePm = (wmean^2)/(4.23^2);
        w = gammaDistr(1:tmax, wmean, shapePm);
    case 5
        % EVD distribution from van Kerkhove 2015
        wmean = 15.3*meanCh; shapePm = (wmean^2)/(9.3^2);
        w = gammaDistr(1:tmax, wmean, shapePm);
    case 6
        % Geometric distribution (prob succ = 1/mean)
        wmean = 20*meanCh; w = geomDistr(1:tmax, wmean);
    case 7
        % Uniform distribution, nWin around mean
        wmean = 10*meanCh; nWin = 7; 
        w = deltaDistr(tmax, wmean, nWin);
    case 8
        % Flare-up from secondary transmission (Lee 2019)
        wmean = 10*meanCh; shapePm = 20;
        w = bimodalGamDistr(1:tmax, wmean, shapePm);
    case 9
        % Uniform distribution looking backwards
        wshift = 20*meanCh; w = zeros(1, tmax); 
        w(1:wshift) = 1/wshift;
end


%% Compute various generation time distributions

%------------------------------------------------------------------
% Gamma distribution with shape and mean
function w = gammaDistr(x, wmean, shapePm)

% Scale parameter based on mean = shapePm*scalePm
scalePm = wmean/shapePm; 
% Gamma probabilities as difference in CDFs
w = gamcdf(x, shapePm, scalePm) - gamcdf(x-1, shapePm, scalePm);

%------------------------------------------------------------------
% Geometric distribution, pSucc is prob success (1/mean)
function w = geomDistr(x, wmean)

% Set pSucc to same dimension as x domain
p1 = 1/wmean; pSucc = p1*ones(size(x));

% Geometric distribution starting at 1 (vs 0)
w = realpow((1-pSucc), (x-1)); w = p1*w;

%------------------------------------------------------------------
% Uniform distribution, window of mass (nWin) on mean
function w = deltaDistr(xmax, wmean, nWin)

% Distribution defined over all time
w = zeros(1, xmax);

% No. time points either side of omega
a = round(wmean) - (nWin - 1)/2;
b = round(wmean) + (nWin - 1)/2;

% Probability over this window
w(a:b) = 1/nWin;

%------------------------------------------------------------------
% Bimodal gamma mixture (fixed shape for both modes)
function w = bimodalGamDistr(x, wmean, shapePm)

% Scale parameter based on mean = shapePm*scalePm
scalePm = wmean/shapePm; ratePm = 1/scalePm;
% Gamma (Erlang) probabilities
w = -log(factorial(shapePm-1)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(x) - ratePm*x;
w = exp(w);

% Flare-up second distribution - larger mean
flareMean = 40; ratePm2 = shapePm/flareMean;
% Log probabilities of second distribution
w2 = -log(factorial(shapePm-1)) + shapePm*log(ratePm2) +...
    (shapePm-1)*log(x) - ratePm2*x; 

% Combined distribution and invert log
w2 = exp(w2); w = w + w2; w = w/sum(w);