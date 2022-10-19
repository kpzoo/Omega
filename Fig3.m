% Compare r, R and omega with changes due to interventions
clearvars; clc; close all; 

% Assumptions and notes
% - compare epidemics with different generation times 
% - ranking changes with interventions

% Directory and where saving
thisDir = cd; saveFol = 'results/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Decide type of generation time distribution (base, pick from set)
epiNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'EVD', 'Geometric',...
    'Uniform', 'Bimodal', 'Window'};
scenNam = {'control', 'square-wave', 'cascade', 'boom-bust', 'filtered',...
    'waves', 'noise valley', 'boom-bust-boom', 'rising', 'falling'};

%% Ebola simulation with two groups and piecewise-constant R

% Times of runs and epidemic type 
epiNo = 5; nday = 201; I0 = 2; wdel = 0.6; 
% Window, generation time and R changes
delta = 30; delR = 0.8; scenNo = 1; 

% True R for both scenarios
simVals1.Rch = [2.2 1.2]; simVals1.tch = 100; 
simVals2.Rch = delR*simVals1.Rch; simVals2.tch = 100; 
simVals = {simVals1, simVals2};

% Main simulation and estimation
[sim1, est1] = variantOmegaSim(epiNo, scenNo, nday, I0, delta, simVals, wdel);


%% COVID simulation with two groups and sinusoidal R

% Times of runs and epidemic type 
epiNo = 4; nday = 201; I0 = 2; wdel = 0.5; 
% Window, generation time and R changes
delta = 15; delR = 0.8; scenNo = 6; 

% True R for both scenarios
simVals1.Rch = [1.2 0.8]; simVals1.tch = 4; 
simVals2.Rch = delR*simVals1.Rch; simVals2.tch = 4; 
simVals = {simVals1, simVals2};

% Main simulation and estimation
[sim2, est2] = variantOmegaSim(epiNo, scenNo, nday, I0, delta, simVals, wdel);

%% Plot estimates on single figure

% EVD (left) and COVID (right)
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
% Axis limits and changes
t1 = [30 sim1.tday(end)]; t2 = [15 sim2.tday(end)];
nw1 = length(sim1.meanPmFrac); nw2 = length(sim2.meanPmFrac);

% Incidence plots
subplot(7, 2, 1);
plot(sim1.tday, sim1.Iday(1, :), 'b', 'LineWidth', 2); hold on;
plot(sim1.tday, sim1.Iday(2, :), 'r', 'LineWidth', 2);
box off; grid off; xlim(t1); hold off;
ylabel('$I_t$', 'FontSize', fnt); title('EVD');
subplot(7, 2, 2); 
plot(sim2.tday, sim2.Iday(1, :), 'b', 'LineWidth', 2); hold on;
plot(sim2.tday, sim2.Iday(2, :), 'r', 'LineWidth', 2);
box off; grid off; xlim(t2); hold off;
ylabel('$I_t$', 'FontSize', fnt); title('COVID-19');

% rt estimates
subplot(7, 2, [3 5]);
stairs(sim1.tday, sim1.rtrue(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t1);
stairs(sim1.tday, sim1.rtrue(2, :), 'k-', 'LineWidth', 2);
plotCIRaw(sim1.tday', est1.rmean(1, :)', est1.rlow(1, :)', est1.rhigh(1, :)', 'b');
plotCIRaw(sim1.tday', est1.rmean(2, :)', est1.rlow(2, :)', est1.rhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{r}_{t}$', 'FontSize', fnt);

subplot(7, 2, [4 6]);
plot(sim2.tday, sim2.rtrue(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t2);
plot(sim2.tday, sim2.rtrue(2, :), 'k-', 'LineWidth', 2);
plotCIRaw(sim2.tday', est2.rmean(1, :)', est2.rlow(1, :)', est2.rhigh(1, :)', 'b');
plotCIRaw(sim2.tday', est2.rmean(2, :)', est2.rlow(2, :)', est2.rhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{r}_{t}$', 'FontSize', fnt);

% Rt estimates
subplot(7, 2, [7 9]);
stairs(sim1.tday, sim1.Rtrue(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t1);
stairs(sim1.tday, sim1.Rtrue(2, :), 'k--', 'LineWidth', 2);
plotCIRaw(sim1.tday', est1.Rmean(1, :)', est1.Rlow(1, :)', est1.Rhigh(1, :)', 'b');
plotCIRaw(sim1.tday', est1.Rmean(2, :)', est1.Rlow(2, :)', est1.Rhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{R}_{t}$', 'FontSize', fnt);

subplot(7, 2, [8 10]);
plot(sim2.tday, sim2.Rtrue(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t2);
plot(sim2.tday, sim2.Rtrue(2, :), 'k-', 'LineWidth', 2);
plotCIRaw(sim2.tday', est2.Rmean(1, :)', est2.Rlow(1, :)', est2.Rhigh(1, :)', 'b');
plotCIRaw(sim2.tday', est2.Rmean(2, :)', est2.Rlow(2, :)', est2.Rhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{R}_{t}$', 'FontSize', fnt);

% Omega estimates
subplot(7, 2, [11 13]);
plot(sim1.tday, sim1.omega(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t1);
plot(sim1.tday, sim1.omega(2, :), 'k-', 'LineWidth', 2);
plotCIRaw(sim1.tday', est1.Ommean(1, :)', est1.Omlow(1, :)', est1.Omhigh(1, :)', 'b');
plotCIRaw(sim1.tday', est1.Ommean(2, :)', est1.Omlow(2, :)', est1.Omhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{\Omega}_{t}$', 'FontSize', fnt);
xlabel('$t$ (d)', 'FontSize', fnt);

subplot(7, 2, [12 14]);
plot(sim2.tday, sim2.omega(1, :), 'k-', 'LineWidth', 2);
hold on; xlim(t2);
plot(sim1.tday, sim2.omega(2, :), 'k-', 'LineWidth', 2);
plotCIRaw(sim2.tday', est2.Ommean(1, :)', est2.Omlow(1, :)', est2.Omhigh(1, :)', 'b');
plotCIRaw(sim2.tday', est2.Ommean(2, :)', est2.Omlow(2, :)', est2.Omhigh(2, :)', 'r');
grid off; box off; hold off;
ylabel('$\hat{\Omega}_{t}$', 'FontSize', fnt);
xlabel('$t$ (d)', 'FontSize', fnt);


