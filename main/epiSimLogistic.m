% Simulate epidemics via renewal model and offsets
function [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimLogistic(scenNo, nday, distNo, simVals)

% Assumptions and notes
% - only works on 2 demes
% - allows for importations between demes
% - presents functions that cancel out
% - removes 'burn-in' of first 20 days, epidemic size < 100
% - various R trajectories adns SI distributions specified

% Possible scenarios available - must match calling function
scenNam = {'step change', 'sine change'};
disp(['True R scenario: ' scenNam{scenNo}]);

% Parameters for changepoints in R trajectory
Rch = simVals.Rch; tch = simVals.tch;
% Variable for true R
Rtrue = zeros(1, nday);

% Functions for scenarios: R on a daily basis
switch(scenNo)
     case 1
        % Either a bottleneck or reverse
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(2)+1:end) = Rch(1);
        if simVals.demeNo == 1
            % Rapid control that recovers
            Rtrue(tch(1)+1:tch(2)) = Rch(1) - Rch(2);
        else
            % Rapid increase then some control
            Rtrue(tch(1)+1:tch(2)) = Rch(1) + Rch(2);
        end
    case 2
        % Either a sine or negative sine
        if simVals.demeNo == 1
            Rtrue = Rch(1) + Rch(2)*sind(tch(1)*(1:nday));
        else
            Rtrue = Rch(1) - Rch(2)*sind(tch(1)*(1:nday));
        end
        
end

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo}; disp(['True SI scenario: ' distChoice]);
% Set serial interval parameters
distvals.type = distNo; 

% Hyerparameters of serial distribution
switch(distNo)
    case 1
        % Geometric distribution - no parameter
        distvals.pm = []; distvals.omega = 15.3;
    case 2
        % Gamma distribution - shape parameter
        %distvals.pm = (1/0.65)^2; distvals.omega = 6.5;
        distvals.pm = 2.7066; distvals.omega = 15.3;
    case 3
        % Delta distribution - odd window around mean
        distvals.pm = 7; distvals.omega = 15.3;
    case 4
        % Two Gamma distributions for flare-up
        distvals.pm = 45; distvals.omega = 15.3;
end

% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Daily incidence and infectiousness
Iday = zeros(1, nday); Lam = Iday; 
% Initialise epidemic and warning
Iday(1) = 10; Iwarn = 0; 

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = sum(Iday(i-1:-1:1).*Pomegat);    
    % Renewal incidence
    Iday(i) = poissrnd(Rtrue(i)*Lam(i));
end

% Remove start-up 20 days
idz = 20:nday; tday = idz;
% Adjusted vectors - including tday
Iday = Iday(idz); Rtrue = Rtrue(idz); Lam = Lam(idz);

% Remove small epidemics
if sum(Iday) < 600
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end
