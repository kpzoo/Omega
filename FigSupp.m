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

% Apply changing omega windows and transformations with canonical windows
lset = 6; dset = round(linspace(10, 50, lset)); 

%% Ebola simulation with piecewise-constant R

% Times of runs and epidemic type 
epiNo = 5; nday = 201; I0 = 2; scenNo = 2;
% Window and generation time changes
delta = 30; wdel = 0.4; wtype = 0;

% Main simulation and estimation
[sim1, est1] = singleOmegaSim(epiNo, scenNo, nday, I0, delta, wdel, wtype);

% Smoothed incidence based growth rate
Iemp1 = smoothdata([0 sim1.Iday]); 
% Smooth and adjust for delay
remp1 = smoothdata(diff(log(Iemp1)), 'sgolay', 30, 'Degree', 2);

% Consider mean estimates under various windows
Omk1 = cell(1, lset); Omd1 = zeros(lset, nday); 
for i = 1:lset
    [Omd1(i, :), Omk1{i}, remp1] = meanWinMetric(sim1.Iday, nday, dset(i), remp1);
end

%% COVID simulation with sinusoidal R

% Times of runs and epidemic type 
epiNo = 4; nday = 201; I0 = 2; scenNo = 6;
% Window and generation time changes
delta = 15; wdel = 0.4; wtype = 1;

% Main simulation and estimation
[sim2, est2] = singleOmegaSim(epiNo, scenNo, nday, I0, delta, wdel, wtype);

% Smoothed incidence based growth rate
Iemp2 = smoothdata([0 sim2.Iday]); 
% Smooth and adjust for delay
remp2 = smoothdata(diff(log(Iemp2)), 'sgolay', 30, 'Degree', 2);

% Consider mean estimates under various windows
Omk2 = cell(1, lset); Omd2 = zeros(lset, nday); 
for i = 1:lset
    [Omd2(i, :), Omk2{i}, remp2] = meanWinMetric(sim2.Iday, nday, dset(i), remp2);
end

%% Compare windows of estimates 

% Axis limits and changes
t1 = [30 sim1.tday(end)]; t2 = [15 sim2.tday(end)];
nw1 = length(sim1.meanPmFrac); nw2 = length(sim2.meanPmFrac);

% Simple plot of true R and estimates
figure('Renderer', 'painters', 'Position', [10 10 800 800]);
subplot(2, 1, 1);
yyaxis('left'); hold on;
plot(sim1.tday, ones(1, nday), '--', 'Color', grey2, 'LineWidth', 0.5);
plot(sim1.tday, sim1.Rtrue, 'k-', 'LineWidth', 2);
plot(sim1.tday, est1.Ommean, 'b-', 'LineWidth', 2);
box off; grid off; xlim(t1); hold off; h = gca;
h.YAxis(1).TickLabelColor = 'k'; h.YAxis(1).Label.Color = 'k';
ylabel('$R, \, \Omega$', 'FontSize', fnt);
yyaxis('right'); hold on;
plot(sim1.tday, zeros(1, nday), '--', 'Color', grey2, 'LineWidth', 0.5);
plot(sim1.tday, remp1, 'r-', 'LineWidth', 2);
sval = range(remp1)*[sim1.wmeans sim1.wmeans(end)]/max(sim1.wmeans);
stairs([sim1.wchtime t1(end)], sval, '--', 'Color', grey2, 'LineWidth', 2);
box off; grid off; xlim(t1); hold off; h = gca;
h.YAxis(2).TickLabelColor = 'k'; h.YAxis(2).Label.Color = 'k';
xlabel(['EVD: $t \, | \, g_0$ = ' num2str(round(sim1.wmeans(1), 1))], 'FontSize', fnt);
ylabel('$r$', 'FontSize', fnt);
subplot(2, 1, 2);
yyaxis('left'); hold on;
plot(sim1.tday, ones(1, nday), '--', 'Color', grey2, 'LineWidth', 0.5);
plot(sim2.tday, sim2.Rtrue, 'k-', 'LineWidth', 2);
plot(sim2.tday, est2.Ommean, 'b-', 'LineWidth', 2);
box off; grid off; xlim(t2); hold off; h = gca;
h.YAxis(1).TickLabelColor = 'k'; h.YAxis(1).Label.Color = 'k';
ylabel('$R, \, \Omega$', 'FontSize', fnt);
yyaxis('right'); hold on;
plot(sim1.tday, zeros(1, nday), '--', 'Color', grey2, 'LineWidth', 0.5);
plot(sim2.tday, remp2, 'r-', 'LineWidth', 2);
sval = range(remp2)*[sim2.wmeans sim2.wmeans(end)]/max(sim2.wmeans);
stairs([sim2.wchtime t2(end)], sval, '--', 'Color', grey2, 'LineWidth', 2);
box off; grid off; xlim(t2); hold off; h = gca;
h.YAxis(2).TickLabelColor = 'k'; h.YAxis(2).Label.Color = 'k';
xlabel(['COVID: $t \, | \, g_0$ = ' num2str(round(sim2.wmeans(1), 1))], 'FontSize', fnt);
ylabel('$r$', 'FontSize', fnt);

% Mean estimate figure for Ebola
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
for i = 1:lset
    subplot(2, 3, i);
    plot(sim1.tday, sim1.Rtrue, 'k-', 'LineWidth', 2);
    hold on;
    plot(sim1.tday, Omk1{i}(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
    plot(sim1.tday, Omk1{i}(1, :), 'Color', 'r', 'LineWidth', 2);
    plot(sim1.tday, Omk1{i}(end, :), 'Color', 'g', 'LineWidth', 2);
    plot(sim1.tday, Omd1(i, :), 'Color', 'b', 'LineWidth', 2);
    plot(sim1.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
    box off; grid off; xlim(t1); hold off; ylim([0 8]);
    ylabel(['$\delta$ = ' num2str(dset(i)) ', ' '$\alpha$ = ' num2str(dset(i)/2)], 'FontSize', fnt);
    if i == 3 || 6
        xlabel(['$t | g_0$ = ' num2str(round(sim1.wmeans(1), 1))], 'FontSize', fnt);
    end
    if i == 1
        legend('$R$ (true)','','','','$R(\alpha)| \sigma$', '$1+r\alpha$',...
            '$e^{r\alpha}$', '$\Omega | \delta$');
    end
end

% Mean estimate figure for COVID
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
for i = 1:lset
    subplot(2, 3, i);
    plot(sim2.tday, sim2.Rtrue, 'k-', 'LineWidth', 2);
    hold on;
    plot(sim2.tday, Omk2{i}(2:end-1, :), 'Color', grey1, 'LineWidth', 2);
    plot(sim2.tday, Omk2{i}(1, :), 'Color', 'r', 'LineWidth', 2);
    plot(sim2.tday, Omk2{i}(end, :), 'Color', 'g', 'LineWidth', 2);
    plot(sim2.tday, Omd2(i, :), 'Color', 'b', 'LineWidth', 2);
    plot(sim2.tday, ones(1, nday), 'k--', 'LineWidth', 1); 
    box off; grid off; xlim(t1); hold off; ylim([0 8]);
    ylabel(['$\delta$ = ' num2str(dset(i)) ', ' '$\alpha$ = ' num2str(dset(i)/2)], 'FontSize', fnt);
    if i == 3 || 6
        xlabel(['$t | g_0$ = ' num2str(round(sim2.wmeans(1), 1))], 'FontSize', fnt);
    end
    if i == 1
        legend('$R$ (true)','','','','$R(\alpha)| \sigma$', '$1+r\alpha$',...
            '$e^{r\alpha}$', '$\Omega | \delta$');
    end
end

%% Plot estimates on single figure

% EVD (left) and COVID (right)
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
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
