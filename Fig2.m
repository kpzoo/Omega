% Deterministic exploration of r, R and omega (Fig 2)
clearvars; clc; close all; 

% Assumptions and notes
% - equations for simple exp growth or decline
% - compare how growth rates and reproduction numbers behave


% Directory and where saving
thisDir = cd; saveFol = 'results/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);


%% Examine how r implies differing R orders but not omega

% Define growth rates r
nr = 200; r = linspace(-0.15, 0.15, nr);
% Infectiousness durations and R
gammas = [12 10 8]; leng = length(gammas);
Rs1 = zeros(leng, nr); Rs2 = Rs1;
for i = 1:leng
    % Exp generation times
    Rs1(i, :) = 1 + r*gammas(i);
    % Delta generation times
    Rs2(i, :) = exp(r*gammas(i));
end

% Compute omega values at del
del = 20; Rs1(Rs1 <= 0) = 0; Rs2(Rs2 <= 0) = 0; 
omega = sqrt((2*r*del)./(1 - exp(-2*r*del)));

% Examine variation with window del
dels = 5:5:30; ndel = length(dels);
omegaD = zeros(ndel, nr);
for i = 1:ndel
    omegaD(i, :) = sqrt((2*r*dels(i))./(1 - exp(-2*r*dels(i))));
end

% Changing R across time
t = 0:30; lent = length(t); tch = median(t);
Rch = 4*ones(1, lent); Rch(tch+1:end) = 3;
% Cause inversion of r and omega via gamma
gamma = 10*ones(1, lent); gamma(tch+1:end) = 5;
rch = (Rch - 1)./gamma; omch = sqrt((2*rch*del)./(1 - exp(-2*rch*del)));

% Publishable figure
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
subplot(6, 2, [1 3 5]);
plot(r, Rs1, 'LineWidth', 2);
hold on;
plot(r, omega, 'k', 'LineWidth', 2);
plot([0 0], [6 1], '--', 'Color', grey2, 'LineWidth', 1);
plot([0 r(1)], [1 1], '--', 'Color', grey2, 'LineWidth', 1);
plot(-[0.1 0.1], [0 6], '--', 'Color', grey1, 'LineWidth', 1);
plot([0.1 0.1], [0 6], '--', 'Color', grey1, 'LineWidth', 1);
hold off; grid off; box off; 
ylim([0 6]); xlim([r(1) r(end)]);
xlabel('$r$ (days$^{-1}$)', 'FontSize', fnt);  ylabel('$R$', 'FontSize', fnt);
leg = legend('', '', '', '$\Omega | \delta = 20$', 'Location', 'best', 'FontSize', fnt); 
leg.Box = 'off'; set(gca,'XAxisLocation','top');

subplot(6, 2, [2 4 6]);
plot(r, Rs2, 'LineWidth', 2);
hold on;
plot(r, omega, 'k', 'LineWidth', 2);
plot([0 0], [6 1], '--', 'Color', grey2, 'LineWidth', 1);
plot([0 r(1)], [1 1], '--', 'Color', grey2, 'LineWidth', 1);
plot(-[0.1 0.1], [0 6], '--', 'Color', grey1, 'LineWidth', 1);
plot([0.1 0.1], [0 6], '--', 'Color', grey1, 'LineWidth', 1);
hold off; grid off; box off; 
ylim([0 6]); xlim([r(1) r(end)]);
xlabel('$r$ (days$^{-1}$)', 'FontSize', fnt);  ylabel('$R$', 'FontSize', fnt);
leg = legend('', '', '', '$\Omega | \delta = 20$', 'Location', 'best', 'FontSize', fnt); 
leg.Box = 'off'; set(gca,'XAxisLocation','top');

subplot(6, 2, [7 9 11]);
plot(r, omegaD(2:end-1, :), 'Color', grey1, 'LineWidth', 2); 
hold on;
plot(r, omegaD(1, :), 'b', 'LineWidth', 2); 
plot(r, omegaD(end, :), 'r', 'LineWidth', 2); 
plot([0 0], [0 1], '--', 'Color', grey2, 'LineWidth', 1);
plot([0 r(1)], [1 1], '--', 'Color', grey2, 'LineWidth', 1);
box off; grid off; 
xlim([r(1) r(end)]); ylim([0 3]);
xlabel('$r$ (days$^{-1}$)', 'FontSize', fnt); 
ylabel('$\Omega | \delta \in [5, 30]$', 'FontSize', fnt);

subplot(6, 2, 8);
stairs(t, Rch, 'b', 'LineWidth', 2);
ylabel('$R$', 'FontSize', fnt);
box off; grid off;
subplot(6, 2, 10);
stairs(t, rch, 'r', 'LineWidth', 2);
ylabel('$r$ (days$^{-1}$)', 'FontSize', fnt);
box off; grid off;
subplot(6, 2, 12);
stairs(t, omch, 'k', 'LineWidth', 2);
box off; grid off;
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$\Omega | \delta = 20$', 'FontSize', fnt);

