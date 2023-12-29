function ts_SC = figures()

% estSC_L2L1_Kendall: [90×90×26 double]
% estSC_L2L1_MSE: [90×90×26 double]
% estSC_L2L1_R: [90×90×26 double]
% kendall_L2L1_ts_SC: [26×4 double]
% maxCorr_L2L1_ts_SC: [26×4 double]
% mse_L2L1_ts_SC: [26×4 double]


group = 'Ctrl';
corrType = 'Pearson';

dataFolder = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC';
copyDir = '/Users/farrasabdelnour/Dropbox/Research/Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);

% Load SC data
mtxSC = zeros(90,90,length(FCdir));
for ii=1:length(FCdir)
    SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(ii).name(10:14) '*sc.csv']);
    tt = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);
    mtxSC(:,:,ii) = tt(1:90,1:90);
end

corr_vs_meanSC = zeros(size(mtxSC,3) , 3);
meanSC = mean( mtxSC , 3);
meanSC = meanSC(:);
%stdSC = std(nonzeros(meanSC(:)));

ts_SC = load([dataFolder filesep 'ts_SC_L1L2_data.mat']);

% Find corr between empirical SC and mean empirical SC
for ii=1:3
    for jj=1:size(mtxSC,3)
        %tt = meanSC;
        %tt = ts_SC.estSC_L2L1_Kendall(:,:,jj);
        tt = mtxSC(:,:,jj);
        tt = tt(:);
        corr_vs_meanSC(jj,1) = corr(meanSC , tt , 'type' , 'Kendall');
        
        %tt = ts_SC.estSC_L2L1_R(:,:,jj);
        tt = mtxSC(:,:,jj);
        tt = tt(:);
        corr_vs_meanSC(jj,2) = corr(meanSC , tt );
        
        %tt = ts_SC.estSC_L2L1_MSE(:,:,jj);
        tt = mtxSC(:,:,jj);
        tt = tt(:);
        corr_vs_meanSC(jj,3) = rmse(meanSC,tt);
    end
end

% Calculate pairwise correlation
for jj = 1:size(mtxSC,3)
    for ii = 1:size(mtxSC,3)

    tt1 = mtxSC(:,:,jj);
    tt2 = mtxSC(:,:,ii);    
    pairwise_kendall(jj,ii) = corr(tt1(:),tt2(:), 'type', 'Kendall');
    pairwise_pearson(jj,ii) = corr(tt1(:),tt2(:));
    pairwise_mse(jj,ii) = rmse(tt1(:),tt2(:));
    end
end


%% Kendall tau figures
figure; 
subplot(2,2,1); 
imagesc(mean(ts_SC.estSC_L2L1_Kendall,3)); colorbar;
title('mean est SC, Kendall');

subplot(2,2,2); 
imagesc(std(ts_SC.estSC_L2L1_Kendall,[],3)); colorbar;
title('std SC, Kendall');

subplot(2,2,3); 
mtx_corr = zeros(size(mtxSC,3) , 5);
mtx_corr(:,1) = ts_SC.kendall_L2L1_ts_SC(:,1);
mtx_corr(:,2) = corr_vs_meanSC(:,1); % Kendall is first column of corr_vs_meanSC
mtx_corr(:,3:end) = ts_SC.kendall_L2L1_ts_SC(:,2:end);
imagesc(mtx_corr); colorbar;
set(gca,'xticklabel',{'Est SC', 'Mean SC','log','L2','L1'} , 'fontsize' , 6);
title('Kendall tau')

subplot(2,2,4); 
clims = [ 0 1];
imagesc(pairwise_kendall , clims ); colorbar
title('Pairwise Kendall tau');

%% Pearson figures
figure; 
subplot(2,2,1); 
imagesc(mean(ts_SC.estSC_L2L1_R,3)); colorbar;
title('mean est SC, Pearson R');

subplot(2,2,2); 
imagesc(std(ts_SC.estSC_L2L1_R,[],3)); colorbar;
title('std SC, Pearson R');

subplot(2,2,3); 
mtx_corr = zeros(size(mtxSC,3) , 5);
mtx_corr(:,1) = ts_SC.maxCorr_L2L1_ts_SC(:,1);
mtx_corr(:,2) = corr_vs_meanSC(:,2); % Pearson is second column of corr_vs_meanSC
mtx_corr(:,3:end) = ts_SC.maxCorr_L2L1_ts_SC(:,2:end);
imagesc(mtx_corr); colorbar;
set(gca,'xticklabel',{'Est SC', 'Mean SC','log','L2','L1'} , 'fontsize' , 6);
title('Pearson R')

subplot(2,2,4); 
clims = [0 1];
imagesc(pairwise_pearson , clims); 
colorbar
title('Pairwise Pearson');

%% MSE figures
figure; 
subplot(2,2,1); 
imagesc(mean(ts_SC.estSC_L2L1_MSE,3)); colorbar;
title('mean est SC, MSE');

subplot(2,2,2); 
imagesc(std(ts_SC.estSC_L2L1_MSE,[],3)); colorbar;
title('std SC, MSE');

subplot(2,2,3); 
mtx_corr = zeros(size(mtxSC,3) , 5);
mtx_corr(:,1) = ts_SC.mse_L2L1_ts_SC(:,1);
mtx_corr(:,2) = corr_vs_meanSC(:,3); 
mtx_corr(:,3:end) = ts_SC.mse_L2L1_ts_SC(:,2:end);
imagesc(mtx_corr); colorbar
set(gca,'xticklabel',{'Est SC', 'Mean SC','log','L2','L1'} , 'fontsize' , 6);
title('MSE')

subplot(2,2,4); 
imagesc(pairwise_mse ); colorbar;
title('Pairwise MSE');




