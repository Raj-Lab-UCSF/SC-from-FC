
numSubj = 26; % Number of NYU healthy subjects
wstdn = 0.001:0.005:0.1; % Noise added to SC intial guess
flagsin.fig = 0;

maxCorr_l1_FC = zeros(numSubj,length(wstdn));
ming_dist_l1_FC = zeros(numSubj,length(wstdn),1);
%Rcoor_l1_FC = zeros(numSubj,length(wstdn),length(wstdn),length(wstdn));
%gdistcoor_l1_FC = zeros(numSubj,2);
Rcoor_l1_FC = {};
gdistcoor_l1_FC = {};

%% L1, 1 - abs(FC_ij) instead of time series input, thresholded FC instead of SC as an initial w_0
tstart = tic;
for kk = 1:numSubj
    for ii=1:length(wstdn)
        [maxCorr_l1_FC(kk,ii),ming_dist_l1_FC(kk,ii), Rcoor , gdistcoor] = runjob_l1_FC( kk, wstdn(ii) , flagsin);
        Rcoor_l1_FC(kk,ii).R = Rcoor;
        gdistcoor_l1_FC(kk,ii).gdist = gdistcoor; 
    end
end
tend = toc(tstart);

%% Do L1/L2 1 - abs(FC_ij) instead of time series input, thresholded FC instead of SC as an initial w_0


%% Do L1 with time series for input, thresholded FC instead of SC for initial w_0

%% Do L1/L2 with time series for input, thresholded FC instead of SC for initial w_0

%% Redo intial work, L1 with time series for input, binarized SC or binarized mean SC as initial w_0

%% Redo intial work, L1/L2 with time series for input, binarized SC or binarized mean SC as initial w_0

%% Run initial work in paper

figure;
subplot(1,2,1)
plot(maxCorr);
title('Max R')
subplot(1,2,2)
plot(ming_dist);
title('Min dist')

tstart = tic;
