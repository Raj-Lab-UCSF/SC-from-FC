% These MM variations use the NYU data
%
% Noise added to binarized SC varies from 0.001 to 0.1 in steps of 0.05

flagsin.fig = 0;
flagsin.save = 1;
saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';

group = 'Ctrl';
corrType = 'Pearson';
thresh = 0.001;
flagsin.SC = 1;
flagsin.FC = 0;

copyDir = '/Users/farrasabdelnour/Dropbox/Research/Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(1).name(10:14) '*sc.csv']);

saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');

%SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);
%SCMtx = SCMtx(1:90,1:90);

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');

% Load SC data
mtxSC = zeros(90,90,length(FCdir));
for ii=1:length(FCdir)
    SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(ii).name(10:14) '*sc.csv']);
    tt = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);
    mtxSC(:,:,ii) = tt(1:90,1:90);
end

meanSC = mean( mtxSC , 3);
stdSC = std(nonzeros(meanSC(:)));
meanSC = meanSC .* (meanSC > stdSC/2);

wstdn = 0; % Noise added to SC intial guess
subj = 1:26; %26; % Number of NYU subjects

%% L2L1, w_0 mean of all SCs, with time series

flagsin.FC = 0;
flagsin.SC = 1;

maxCorr_L2L1_ts_SC = zeros( length(subj) , 4 );
%ming_distL2L1_ts_SC = zeros( length(subj) , 4 );
kendall_L2L1_ts_SC = zeros( length(subj) , 4 );
mse_L2L1_ts_SC = zeros( length(subj) , 4 );

estSC_L2L1_R = zeros(90,90,length(subj));
estSC_L2L1_Kendall = zeros(90,90,length(subj));
estSC_L2L1_MSE = zeros(90,90,length(subj));

tStart = tic;
for ii = 1:length(subj)
%    for jj = 1:length(wstdn)
     %tic
     [ maxCorr_L2L1_ts_SC(ii,:) , kendall_L2L1_ts_SC(ii,:) , mse_L2L1_ts_SC(ii,:) , estSC_L2L1_R(:,:,ii) , estSC_L2L1_Kendall(:,:,ii) , estSC_L2L1_MSE(:,:,ii)] = runjob_l1l2( ii , mtxSC(:,:,ii) , meanSC , flagsin);
     %toc
%    end
end
tEnd = toc(tStart);

disp(tEnd/60/60)

save ts_SC_L1L2_data.mat maxCorr_L2L1_ts_SC kendall_L2L1_ts_SC mse_L2L1_ts_SC estSC_L2L1_R estSC_L2L1_Kendall estSC_L2L1_MSE;

% figure;
% imagesc( maxCorrL2L1_ts_SC )
% title('R L2L1 with time series and binarized SC') 
% saveas(gcf , [saveFigs filesep 'maxCorrL2L1_ts_SC'] , 'jpg');
% 
% figure;
% imagesc( ming_distL2L1_ts_SC )
% title('Riemann, L2L1 with time series and binarized SC')
% saveas(gcf , [saveFigs filesep 'ming_distL2L1_ts_SC'] , 'jpg');

%% L2L1, w_0 is FC
%% L2L1, w_0 mean of all SCs, with time series

flagsin.fig = 0;
flagsin.FC = 1;
flagsin.SC = 0;

maxCorr_L2L1_ts_FC = zeros( length(subj) , 4 );
kendall_L2L1_ts_FC = zeros( length(subj) , 4 );
mse_L2L1_ts_FC = zeros( length(subj) , 4 );

estSC_L2L1_ts_FC_R = zeros(90,90,length(subj));
estSC_L2L1_ts_FC_Kendall = zeros(90,90,length(subj));
estSC_L2L1_ts_FC_MSE = zeros(90,90,length(subj));

tStart = tic;
for ii=1:length(subj)
    [maxCorr_L2L1_ts_FC(ii,:) , kendall_L2L1_ts_FC(ii,:) , mse_L2L1_ts_FC(ii,:) , estSC_L2L1_ts_FC_R(:,:,ii) , estSC_L2L1_ts_FC_Kendall(:,:,ii) , estSC_L2L1_ts_FC_MSE(:,:,ii)] = runjob_l1l2_ts_FC( ii , mtxSC(:,:,ii) , flagsin);
end
tEnd = toc(tStart);

disp(tEnd/60/60)
save ts_FC_L1L2_data.mat maxCorr_L2L1_ts_FC kendall_L2L1_ts_FC mse_L2L1_ts_FC estSC_L2L1_ts_FC_R estSC_L2L1_ts_FC_Kendall estSC_L2L1_ts_FC_MSE;


