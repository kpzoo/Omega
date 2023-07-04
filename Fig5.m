% Assess prediction error on examples of Fig 3 to get Fig 5
clearvars; clc; close all; tic;

% Assumptions and notes
% - only filtered estimates are used
% - estimate R and Omega and use one-step-ahead predictions
% - consider changes in generation time as done in Fig 3
% - simulate M trajectories and obtain APE and PMSE

% Directory and where saving
thisDir = cd; saveFol = 'batch/'; 
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

% Define number of simulated replicates
M = 200; disp(['Replicate number = ' num2str(M)]);

%% Ebola simulation with piecewise-constant R

% Times of runs and epidemic type 
epiNo = 5; nday = 201; I0 = 2; scenNo = 2;
% Window and generation time changes
deltaEVD = 30; wdel = 0.4; wtype = 0;

% Variables for storing loop outputs
apeREVD = zeros(1, M); apeR0EVD = apeREVD; apeOmEVD = apeREVD; 
pmseREVD = apeREVD; pmseR0EVD = apeREVD; pmseOmEVD = apeREVD;
% Single simulation outputs
Ievd = cell(1, 1); tevd = Ievd; Revd = Ievd; Omevd = Revd; 
R0evd = Revd; IRevd = Ievd; IOmevd = Ievd; IR0evd = Ievd;
IRevdci = Ievd; IOmevdci = Ievd; IR0evdci = Ievd;

parfor i = 1:M
    % Main simulation and estimation
    [sim1, est1] = singleOmegaSimPred(epiNo, scenNo, nday, I0, deltaEVD, wdel, wtype);

    % One-step-ahead scores
    apeREVD(i) = est1.apeR; apeR0EVD(i) = est1.apeR0; apeOmEVD(i) = est1.apeOm;
    pmseREVD(i) = est1.pmseR; pmseR0EVD(i) = est1.pmseR0; pmseOmEVD(i) = est1.pmseOm;

    % Take one example
    if i == 1
        Ievd{i} = sim1.Iday; tevd{i} = sim1.tday;
        Revd{i} = est1.Rmean; Omevd{i} = est1.Ommean; R0evd{i} = est1.Rmean0;
        IRevd{i} = est1.IRpred; IOmevd{i} = est1.IOmpred; IR0evd{i} = est1.IRpred0;
        IRevdci{i} = est1.IRci; IOmevdci{i} = est1.IOmci; IR0evdci{i} = est1.IRci0; 
    end

    disp(['Completed ' num2str(i) ' of ' num2str(M)]);
end


% Timing and data saving
tsim = toc/60; disp(['Run time = ' num2str(tsim)]); tstamp = datetime;
% Saving data and figs name
namstr = ['batch_' num2str(epiNo) '_' num2str(scenNo) '_' num2str(deltaEVD) '_' num2str(M)];

% Remove unnecessary variables if saving
if saveTrue
    cd(saveFol);
    save(['EVD_' namstr '.mat']);
    cd(thisDir);
end


%% COVID simulation with sinusoidal R

% Times of runs and epidemic type 
epiNo = 4; nday = 201; I0 = 2; scenNo = 6;
% Window and generation time changes
deltaCOVID = 15; wdel = 0.6; wtype = 1;

% Variables for storing loop outputs
apeRCOVID = zeros(1, M); apeR0COVID = apeRCOVID; apeOmCOVID = apeRCOVID; 
pmseRCOVID = apeRCOVID; pmseR0COVID = apeRCOVID; pmseOmCOVID = apeRCOVID;

% Single simulation outputs
Icov = cell(1, 1); tcov = Icov; Rcov = Icov; Omcov = Icov; 
R0cov = Icov; IRcov = Icov; IOmcov = Icov; IR0cov = Icov;
IRcovci = Icov; IOmcovci = Icov; IR0covci = Icov;

parfor i = 1:M
    % Main simulation and estimation
    [sim2, est2] = singleOmegaSimPred(epiNo, scenNo, nday, I0, deltaCOVID, wdel, wtype);

    % One-step-ahead scores
    apeRCOVID(i) = est2.apeR; apeR0COVID(i) = est2.apeR0; apeOmCOVID(i) = est2.apeOm;
    pmseRCOVID(i) = est2.pmseR; pmseR0COVID(i) = est2.pmseR0; pmseOmCOVID(i) = est2.pmseOm;

    % Take one example
    if i == 1
        Icov{i} = sim2.Iday; tcov{i} = sim2.tday;
        Rcov{i} = est2.Rmean; Omcov{i} = est2.Ommean; R0cov{i} = est2.Rmean0;
        IRcov{i} = est2.IRpred; IOmcov{i} = est2.IOmpred; IR0cov{i} = est2.IRpred0;
        IRcovci{i} = est2.IRci; IOmcovci{i} = est2.IOmci; IR0covci{i} = est2.IRci0; 
    end

    disp(['Completed ' num2str(i) ' of ' num2str(M)]);
end


% Timing and data saving
tsim = toc/60; disp(['Run time = ' num2str(tsim)]); tstamp = datetime;
% Saving data and figs name
namstr = ['batch_' num2str(epiNo) '_' num2str(scenNo) '_' num2str(deltaCOVID) '_' num2str(M)];

% Remove unnecessary variables if saving
if saveTrue
    cd(saveFol);
    save(['COVID_' namstr '.mat']);
    cd(thisDir);
end

%% Publishable figure summarising prediction results

% Axis limits for EVD and COVID
t1 = [30 tevd{1}(end)]; t2 = [15 tcov{1}(end)]; t0 = 2:nday;

% EVD (left) and COVID (right)
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);

% Misspecified R based predictions
subplot(3, 2, 1);
plot(tevd{1}, Ievd{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IRevd{1}', IRevdci{1}(:, 1), IRevdci{1}(:, 2), 'r'); 
ylabel('$I_t | \hat{R}_t$', 'FontSize', fnt); title('EVD');
box off; grid off; hold off;
subplot(3, 2, 2);
plot(tcov{1}, Icov{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IRcov{1}', IRcovci{1}(:, 1), IRcovci{1}(:, 2), 'r'); 
ylabel('$I_t | \hat{R}_t$', 'FontSize', fnt); title('COVID');
box off; grid off; hold off;

% Omega predictions
subplot(3, 2, 3);
plot(tevd{1}, Ievd{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IOmevd{1}', IOmevdci{1}(:, 1), IOmevdci{1}(:, 2), 'b'); 
ylabel('$I_t | \hat{\Omega}_t$', 'FontSize', fnt); 
box off; grid off; hold off; xlabel('$t$ (d)', 'FontSize', fnt);
subplot(3, 2, 4);
plot(tcov{1}, Icov{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IOmcov{1}', IOmcovci{1}(:, 1), IOmcovci{1}(:, 2), 'b'); 
ylabel('$I_t | \hat{\Omega}_t$', 'FontSize', fnt); 
box off; grid off; hold off; xlabel('$t$ (d)', 'FontSize', fnt);

% Boxplot of prediction quality
subplot(3, 2, 5);
[fR0, x0] = ksdensity(pmseR0EVD);
fR = ksdensity(pmseREVD, x0); fOm = ksdensity(pmseOmEVD, x0);
plot(x0, fR0, 'k', x0, fR, 'r', x0, fOm, 'b', 'LineWidth', 2);
[fR0, x1] = ksdensity(apeR0EVD); hold on;
fR = ksdensity(apeREVD, x1); fOm = ksdensity(apeOmEVD, x1);
plot(x1, fR0, 'k--', x1, fR, 'r--', x1, fOm, 'b--', 'LineWidth', 2);
xlabel('$D(I_t | I_1^{t-1})$', 'FontSize', fnt);
box off; grid off; hold off; ylabel('P($D$)', 'FontSize', fnt);

subplot(3, 2, 6);
[fR0, x0] = ksdensity(pmseR0COVID);
fR = ksdensity(pmseRCOVID, x0); fOm = ksdensity(pmseOmCOVID, x0);
plot(x0, fR0, 'k', x0, fR, 'r', x0, fOm, 'b', 'LineWidth', 2);
[fR0, x1] = ksdensity(apeR0COVID); hold on;
fR = ksdensity(apeRCOVID, x1); fOm = ksdensity(apeOmCOVID, x1);
plot(x1, fR0, 'k--', x1, fR, 'r--', x1, fOm, 'b--', 'LineWidth', 2);
xlabel('$D(I_t | I_1^{t-1})$', 'FontSize', fnt);
box off; grid off; hold off; ylabel('P($D$)', 'FontSize', fnt);



%% Version of figure with histograms


% EVD (left) and COVID (right)
figure('Renderer', 'painters', 'Position', [10 10 1200 1000]);

% Misspecified R based predictions
subplot(6, 2, [1 3]);
plot(tevd{1}, Ievd{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IRevd{1}', IRevdci{1}(:, 1), IRevdci{1}(:, 2), 'r'); 
ylabel('$I_t | \hat{R}_t$', 'FontSize', fnt); title('EVD');
box off; grid off; hold off;
subplot(6, 2, [2 4]);
plot(tcov{1}, Icov{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IRcov{1}', IRcovci{1}(:, 1), IRcovci{1}(:, 2), 'r'); 
ylabel('$I_t | \hat{R}_t$', 'FontSize', fnt); title('COVID');
box off; grid off; hold off;

% Omega predictions
subplot(6, 2, [5 7]);
plot(tevd{1}, Ievd{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IOmevd{1}', IOmevdci{1}(:, 1), IOmevdci{1}(:, 2), 'b'); 
ylabel('$I_t | \hat{\Omega}_t$', 'FontSize', fnt); 
box off; grid off; hold off; xlabel('$t$ (d)', 'FontSize', fnt);
subplot(6, 2, [6 8]);
plot(tcov{1}, Icov{1}, '.', 'Color', 'k', 'MarkerSize', 20);
hold on; xlim(t1);
plotCIRaw(t0', IOmcov{1}', IOmcovci{1}(:, 1), IOmcovci{1}(:, 2), 'b'); 
ylabel('$I_t | \hat{\Omega}_t$', 'FontSize', fnt); 
box off; grid off; hold off; xlabel('$t$ (d)', 'FontSize', fnt);

% Histograms of prediction quality
subplot(6, 2, 9); nbin = 20;
h = histogram(apeR0EVD, nbin, 'Normalization', 'probability');
h.EdgeAlpha = 0.2; h.FaceColor = "none"; h.EdgeColor = 'k';
hold on; box off; grid off;
h = histogram(apeREVD, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'r';
h = histogram(apeOmEVD, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'b';
xlabel('APE', 'FontSize', fnt);

subplot(6, 2, 11); nbin = 20;
h = histogram(pmseR0EVD, nbin, 'Normalization', 'probability');
h.EdgeAlpha = 0.2; h.FaceColor = "none"; h.EdgeColor = 'k';
hold on; box off; grid off;
h = histogram(pmseREVD, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'r';
h = histogram(pmseOmEVD, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'b';
xlabel('PMSE', 'FontSize', fnt);

subplot(6, 2, 10); nbin = 20;
h = histogram(apeR0COVID, nbin, 'Normalization', 'probability');
h.EdgeAlpha = 0.2; h.FaceColor = "none"; h.EdgeColor = 'k';
hold on; box off; grid off;
h = histogram(apeRCOVID, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.FaceColor = "none"; h.FaceColor = 'r';
h = histogram(apeOmCOVID, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'b';
xlabel('APE', 'FontSize', fnt);

subplot(6, 2, 12); nbin = 20;
h = histogram(pmseR0COVID, nbin, 'Normalization', 'probability');
h.EdgeAlpha = 0.2; h.FaceColor = "none"; h.EdgeColor = 'k';
hold on; box off; grid off;
h = histogram(pmseRCOVID, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'r';
h = histogram(pmseOmCOVID, nbin, 'Normalization', 'probability');
h.FaceAlpha = 0.12; h.EdgeColor = "none"; h.FaceColor = 'b';
xlabel('PMSE', 'FontSize', fnt);


