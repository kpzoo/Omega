% Estimate R and omega in the face of varying generation times
clearvars; clc; close all; 

% Assumptions and notes
% - examine true values and estimates (EpiFilter)
% - change R and generation time distribution

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


%% Ebola simulation with piecewise-constant R

% Times of runs and epidemic type 
epiNo = 5; nday = 201; I0 = 2; scenNo = 2;
% Window and generation time changes
delta = 30; wdel = 0.4; wtype = 0;

% Main simulation and estimation
[sim1, est1] = singleOmegaSim(epiNo, scenNo, nday, I0, delta, wdel, wtype);


%% COVID simulation with sinusoidal R

% Times of runs and epidemic type 
epiNo = 4; nday = 201; I0 = 2; scenNo = 6;
% Window and generation time changes
delta = 15; wdel = 0.4; wtype = 1;

% Main simulation and estimation
[sim2, est2] = singleOmegaSim(epiNo, scenNo, nday, I0, delta, wdel, wtype);

%% Plot estimates on single figure

% EVD (left) and COVID (right)
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
% Axis limits and changes
t1 = [30 sim1.tday(end)]; t2 = [15 sim2.tday(end)];
nw1 = length(sim1.meanPmFrac); nw2 = length(sim2.meanPmFrac);

% Incidence plots
subplot(5, 2, 1);
yyaxis('left');
plot(sim1.tday, sim1.Iday, 'k', 'LineWidth', 2); hold on;
h = gca; yval = [0 h.YLim(2)]; 
plot(repmat(sim1.wchtime(2:end), [2 1]), repmat(yval', [1 nw1-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
box off; grid off; xlim(t1); hold off;
h.YAxis(1).TickLabelColor = 'k'; h.YAxis(1).Label.Color = 'k';
ylabel('$I_t$', 'FontSize', fnt); title('EVD');
yyaxis('right');
stairs([sim1.wchtime t1(end)], [sim1.wmeans sim1.wmeans(end)], 'Color', grey2,...
    'LineWidth', 1);
h = gca; h.YAxis(2).TickLabelColor = 'k'; box off; grid off; 

subplot(5, 2, 2); 
plot(sim2.tday, sim2.Iday, 'k', 'LineWidth', 2); hold on;
h = gca; yval = [0 h.YLim(2)]; 
plot(repmat(sim2.wchtime(2:end), [2 1]), repmat(yval', [1 nw2-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
box off; grid off; xlim(t2); hold off;
h.YAxis(1).TickLabelColor = 'k'; h.YAxis(1).Label.Color = 'k';
ylabel('$I_t$', 'FontSize', fnt); title('COVID-19');
yyaxis('right');
stairs([sim2.wchtime t2(end)], [sim2.wmeans sim2.wmeans(end)], 'Color', grey2,...
    'LineWidth', 1);
h = gca; h.YAxis(2).TickLabelColor = 'k'; box off; grid off; 

% Rt estimates
subplot(5, 2, [3 5]);
plot(sim1.tday, sim1.Rtrue, 'k-', 'LineWidth', 2);
hold on; xlim(t1);
plotCIRaw(sim1.tday', est1.Rmean', est1.Rlow', est1.Rhigh', 'r');
h = gca; yval = [0 h.YLim(2)];
plot(repmat(sim1.wchtime(2:end), [2 1]), repmat(yval', [1 nw1-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
plot(sim1.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
grid off; box off; hold off;
ylabel('$\hat{R}_{t}$', 'FontSize', fnt);

subplot(5, 2, [4 6]);
plot(sim2.tday, sim2.Rtrue, 'k-', 'LineWidth', 2);
hold on; xlim(t2);
plotCIRaw(sim2.tday', est2.Rmean', est2.Rlow', est2.Rhigh', 'r');
h = gca; yval = [0 h.YLim(2)];
plot(repmat(sim2.wchtime(2:end), [2 1]), repmat(yval', [1 nw2-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
plot(sim2.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
grid off; box off; hold off;
ylabel('$\hat{R}_{t}$', 'FontSize', fnt);

% Omega estimates
subplot(5, 2, [7 9]);
plot(sim1.tday, sim1.omega, 'k-', 'LineWidth', 2);
hold on; xlim(t1);
plotCIRaw(sim1.tday', est1.Ommean', est1.Omlow', est1.Omhigh', 'b');
h = gca; yval = [0 h.YLim(2)];
plot(repmat(sim1.wchtime(2:end), [2 1]), repmat(yval', [1 nw1-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
plot(sim1.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
grid off; box off; hold off;
ylabel('$\hat{\Omega}_{t}$', 'FontSize', fnt);
xlabel('$t$ (d)', 'FontSize', fnt);

subplot(5, 2, [8 10]);
plot(sim2.tday, sim2.omega, 'k-', 'LineWidth', 2);
hold on; xlim(t2);
plotCIRaw(sim2.tday', est2.Ommean', est2.Omlow', est2.Omhigh', 'b');
h = gca; yval = [0 h.YLim(2)];
plot(repmat(sim2.wchtime(2:end), [2 1]), repmat(yval', [1 nw2-1]), '--',...
    'Color', grey1, 'LineWidth', 1);
plot(sim2.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
grid off; box off; hold off;
ylabel('$\hat{\Omega}_{t}$', 'FontSize', fnt);
xlabel('$t$ (d)', 'FontSize', fnt);

% % Generation time distributions
% subplot(6, 2, 11);
% plot(sim1.tday, sim1.wch, 'LineWidth', 2);
% xlim([0 40]); xlabel('$u$ (d)', 'FontSize', fnt);
% ylabel('$w_u$', 'FontSize', fnt);
% grid off; box off;
% 
% subplot(6, 2, 12);
% plot(sim2.tday, sim2.wch, 'LineWidth', 2);
% xlim([0 40]); xlabel('$u$ (d)', 'FontSize', fnt);
% ylabel('$w_u$', 'FontSize', fnt);
% grid off; box off;
