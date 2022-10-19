% Setup true R scenarios
function Rtrue = scenarioReprod(scenNo, nday, simVals)

% Assumptions and notes
% - various R trajectories for use in simulation

% Possible scenarios available - must match calling function
scenNam = {'control', 'square-wave', 'cascade', 'boom-bust', 'filtered',...
    'waves', 'noise valley', 'boom-bust-boom', 'rising', 'falling'};
%disp(['True R scenario: ' scenNam{scenNo}]);

%% True R number scenarios to simulate

% Parameters for changepoints in R trajectory
Rch = simVals.Rch; tch = simVals.tch;
% Variable for true R
Rtrue = zeros(1, nday);

% Functions for scenarios: R on a daily basis
switch(scenNo)
     case 1
        % Rapidly controlled epidemic
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
    case 2
        % Rapid control that recovers
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
    case 3
        % Three stage control with fluctuations
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:tch(3)) = Rch(3);
        Rtrue(tch(3)+1:end) = Rch(4);
        Rtrue = Rtrue + 0.3*cosd(2*(1:nday));
        Rtrue(Rtrue <= 0) = 0.1;
    case 4
        % Exponential rise and fall
        trise = 1:tch; tfall = tch+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.03*(1:tch)); Rmax = Rtrue(tch);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.015*(tfall - tch));
    case 5
        % Two stage control with filtered noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
    case 6
        % Second (sine) wave dynamics
        Rtrue = Rch(1) + Rch(2)*sind(tch(1)*(1:nday));
    case 7
        % Long period of low R between transmission and noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
        % Check for noise
        if any(Rtrue < 0)
            Rtrue(Rtrue < 0) = 0;
        end
    case 8
        % Exponential rise and fall then rise
        tch1 = tch(1); tch2 = tch(2);
        trise = 1:tch1; tfall = tch1+1:tch2; 
        triseSec = tch2+1:nday; % second wave
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.03*(1:tch1)); Rmax = Rtrue(tch1);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.015*(tfall - tch1));
        % Second wave
        Rtrue(triseSec) = Rtrue(tch2)*exp(0.02*(triseSec - tch2));
    case 9
        % Simple increase in transmissibility all at R > 1
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
    case 10
        % Simple decrease in transmissibility 
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
end