% Empirical example using Hong Kong data (also in R)
clearvars; clc; close all; tic;

% Assumptions and notes
% - estimate R and Omega with updated and wrong serial intervals
% - also provide estimates using growth rates

% Directory and where saving (and loading)
thisDir = cd; saveFol = 'hkdata/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Load incidence data and serial intervals
cd('hkdata');
Idata = readtable('cases_by_infectee_onset_date.csv');
sidata = readtable('serial_intervals.csv');
cd(thisDir);

%% Prepare data from Ali et al 2020

% Format data into time series
Iday = Idata.n'; nday = length(Iday);
tdate = Idata.infectee_onsetDate;
xdate = datenum(tdate); tday = 1:nday;
tdate = datetime(tdate, 'InputFormat', 'yyyy/MM/dd');

% Change times of serial intervals
sichange = unique(sidata.date_set); nsi = length(sichange);
% Data informing different serial intervals
nsi_val = length(sidata.date_set); idsi = zeros(1, nsi_val);
for i = 1:nsi_val
    idsi(i) = find(strcmp(sidata.date_set(i), sichange));
end

% Manually find which span whole period and which are segments
siall = sichange(3); siseg = sichange([4 1 2]);
% Collect serial interval data for whole period
sidata_all = sidata.si(idsi == 3); idch = zeros(1, 3);

% Define ending change times of segments
tch = {'2020/01/22', '2020/01/29', '2020/02/13'};
tch = datetime(tch, 'InputFormat', 'yyyy/MM/dd');
% Collect serial data for each segment in order (manual)
sidata_seg = cell(1, 3); idseg = [4 1 2];
mean_seg = zeros(1, 3); sd_seg = mean_seg;
for i = 1:3
   sidata_seg{i} = sidata.si(idsi == idseg(i));
   % Fit normal distribution
   [mu, sig] = normfit(sidata_seg{i});
   mean_seg(i) = mu; sd_seg(i) = sig; 
   % Change time for serial intervals
   idch(i) = find(tdate == tch(i));
end

% Direct values from Ali et al 2020
meanCh = [7.1, 5.2, 3.0]; sdCh = [5.3, 4.7, 4.1];
% Truncate to last value for SI
nday = idch(end); Iday = Iday(1:nday); tday = tday(1:nday);
tdate = tdate(1:nday); xdate = xdate(1:nday);

%% Estimate R using fixed and varying SI and omega

% Window and changes
delta = 15; nch = length(idch); 

% Main estimation function for R and omega (includes growth rate)
%[sim1, est1] = hkOmega(nday, Iday, delta, meanCh, sdCh.^2, idch);
[sim1, est1] = hkOmega(nday, Iday, delta, mean_seg, sd_seg.^2, idch);

% Timing and data saving
tsim = toc/60; disp(['Run time = ' num2str(tsim)]); tstamp = datetime;
% Saving data and figs name
namstr = ['HK_' num2str(nch) '_' num2str(delta)];

% Examine simple fits to histogram
x = 0:0.01:25; 
figure;
for i = 1:3
    % Gamma parameters
    scalePm = (sd_seg(i)^2)/mean_seg(i); shapePm = mean_seg(i)/scalePm;
    subplot(3, 1, i);
    histogram(sidata_seg{i}, 20, 'Normalization', 'probability');
    hold on;
    plot(x, gampdf(x, shapePm, scalePm), 'r', 'LineWidth', 2);
    hold off; box off; grid off;
end

% Smoothed incidence based growth rate
Iemp = smoothdata([0 Iday]); 
% Smooth and adjust for delay
remp = diff(log(Iemp)); remp = smoothdata(remp);


% Remove unnecessary variables if saving
if saveTrue
    cd(saveFol);
    save([namstr '.mat']);
    cd(thisDir);
end

%% Publishable figure comparing estimates

% Axes for control
ax = zeros(1, 4); xdate = xdate';

figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
for i = 1:3
    % Plot against fixed serial intervals
    ax(i) = subplot(3, 2, 2*i-1);
    plotCIRaw(xdate', est1.Rfix(i, :)', est1.Rlfix(i, :)', est1.Rhfix(i, :)', 'g');
    hold on;
    plotCIRaw(xdate', est1.Ommean', est1.Omlow', est1.Omhigh', 'b');
    plot(xdate, ones(1, nday), 'k--', 'LineWidth', 1); 
    h = gca; 
    plot(xdate(idch(1))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
    plot(xdate(idch(2))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
    hold off; box off; grid off;
    xlim([xdate(1) xdate(end)]);
    datetick('x','dd-mm', 'keepticks');
    ylabel(['$\hat{R}_t$[' num2str(i) '], $\hat{\Omega}_t$'], 'FontSize', fnt);
    if i == 3
        xlabel('$t$ (d), $\delta = 15d$', 'FontSize', fnt);
    end
end

% Time varying serial intervals
ax(4) = subplot(3, 2, 4);
plotCIRaw(xdate', est1.Rmean', est1.Rlow', est1.Rhigh', 'r');
hold on;
plotCIRaw(xdate', est1.Ommean', est1.Omlow', est1.Omhigh', 'b');
plot(xdate, ones(1, nday), 'k--', 'LineWidth', 1); 
h = gca; 
plot(xdate(idch(1))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
plot(xdate(idch(2))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
hold off; box off; grid off;
xlim([xdate(1) xdate(end)]);
datetick('x','dd-mm', 'keepticks');
%linkaxes(ax, 'xy');
ylabel('$\hat{R}_t, \hat{\Omega}_t$', 'FontSize', fnt);

% Incidence and denominators
subplot(3, 2, 2);
yyaxis('left');
stairs(xdate, Iday, 'k', 'LineWidth', 2);
hold on;
plot(xdate, sim1.Lday, 'r-', 'LineWidth', 2);
plot(xdate, sim1. Ipmean, 'b-', 'LineWidth', 2);
h = gca; h.YAxis(1).TickLabelColor = 'k'; h.YAxis(1).Label.Color = 'k';
plot(xdate(idch(1))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
plot(xdate(idch(2))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
hold off; box off; grid off;
xlim([xdate(1) xdate(end)]);
legend('', '$\Lambda_t$', '$M_t$','', '', 'Location','best');
datetick('x','dd-mm', 'keepticks');
ylabel('$I_t, \Lambda_t, M_t$', 'FontSize', fnt);
yyaxis('right');
stairs(xdate([1 idch]), [mean_seg mean_seg(end)], 'Color', grey2, 'LineWidth', 2);
h = gca; h.YAxis(2).TickLabelColor = 'k'; h.YAxis(2).Label.Color = 'k';

% Growth rates
ax(4) = subplot(3, 2, 6);
plotCIRaw(xdate', est1.rmean', est1.rlow', est1.rhigh', 'r');
hold on;
plotCIRaw(xdate', est1.rommean', est1.romlow', est1.romhigh', 'b');
plot(xdate, remp, 'k', 'LineWidth', 2);
plot(xdate, zeros(1, nday), 'k--', 'LineWidth', 1); 
h = gca; 
plot(xdate(idch(1))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
plot(xdate(idch(2))*ones(1, 2), h.YLim, '--', 'Color', grey2, 'LineWidth', 1);
hold off; box off; grid off;
xlim([xdate(1) xdate(end)]);
datetick('x','dd-mm', 'keepticks');
ylabel('$\hat{r}_t | \hat{R}_t, \hat{r}_t | \hat{\Omega}_t$', 'FontSize', fnt);
xlabel('$t$ (d), $\delta = 15d$', 'FontSize', fnt);

% figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
% % Plot against fixed serial intervals
% subplot(2, 2, 1); hold on;
% for i = 1:3    
%     plotCIRaw(tday', est1.Rfix(i, :)', est1.Rlfix(i, :)', est1.Rhfix(i, :)', 'g');
% end
% plotCIRaw(tday', est1.Rmean', est1.Rlow', est1.Rhigh', 'r');
% subplot(2, 2, 3); hold on;
% for i = 1:3    
%     plotCIRaw(tday', est1.Rfix(i, :)', est1.Rlfix(i, :)', est1.Rhfix(i, :)', 'g');
% end
% plotCIRaw(tday', est1.Ommean', est1.Omlow', est1.Omhigh', 'b');
