function [maxCorrArray,kendall_valArray,mse_valArray,w2d,w_kendall,w_mse] = runjob_l1l2_ts_FC(subj , subjSC , flagsin)

% Here we perform L1/L2 penalty for estimating SC
%
% Subject's SC is binarized with added noise to the zero elements of SC
%
% Majorization-minimization algorithm operates on *time series* and *L1/L2*
% on *binarized SC*

group = 'Ctrl';
%corrType = 'Pearson';
thresh = 0.001;
flagsin.SC = 0;
flagsin.FC = 1;

%wstdn = 2^(-3);

copyDir = '/Users/farrasabdelnour/Dropbox/Research//Projects'; 
%FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
%SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);
FCdir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*Cov_FC.mat');
SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);

saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';
%dataLoc = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/Data/MEG/meg_individual_tseries_reordered.nc';
%finfo = ncinfo(dataLoc);
%vardata = ncread( dataLoc , values );

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');
%SCMtx = csvread([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name],1,0);
%SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);
%SCMtx = SCMtx(1:90,1:90);


tsl = load([tsFCDir(subj).folder filesep tsFCDir(subj).name]);

%% Vectorize FC upper triangular
% FC = abs(FC);
% 
fsl = load([FCdir(subj).folder filesep FCdir(subj).name]);
FC = fsl.FCPearson;

absFC = abs(FC);
stdFC = std(nonzeros(absFC(:)));
threshFC = absFC .* (absFC > 2*stdFC);
threshFC = threshFC.';
m  = (1:size(threshFC,1)).' > (1:size(threshFC,2));
w_f  = threshFC(m);
clear m absFC stdFC threshFC;

%% Vectorize empirical SC upper triangle
subjSC = subjSC.';
m  = (1:size(subjSC,1)).' > (1:size(subjSC,2));
w  = subjSC(m);
clear subjSC;

% Collect fMRI time series
mtxFC = zeros(90,190);
for ii=1:90
    tt = eval(['tsl.AAL' num2str(ii,'%2.2d') 'TC']);
    mtxFC(ii,:) =  tt';
end

nSamples = size(mtxFC , 2);

f = 0.001;
a = (f:0.05:2);
b = (f:0.05:2);
c = (f:0.05:2);

% Initialize metrics
maxCorr = 0;
kendall_val = 0;
%ming_dist = inf;
mse_val = inf;

corrTable = zeros(length(a) , length(b) , length(c)); % Pearson corr
%g_dist = zeros(length(a) , length(b) , length(c)); % Riemann 
mse_array = zeros(size(corrTable)); % MSE
kendall_array = zeros(size(corrTable)); % Kendall tau

for ii = 1:length(a)
    %disp(a(ii))
    for jj = 1:length(b)
        %disp(b(jj))
        for kk = 1:length(c)
            %disp(c(kk))
            [wk, stat] = MM_gl_l2l1_ts_FC(mtxFC, a(ii), b(jj),c(kk),w_f,thresh,nSamples);
            wk = wk .* (wk > std(wk)/2); % May remove later
            corrTable(ii,jj,kk) = corr(w,wk);
            kendall_array(ii,jj,kk) = corr(w,wk,'type','Kendall');
            %g_dist(ii,jj,kk) = riemann_dist( squareform(w) , squareform(wk) , flagsin );
            mse_array(ii,jj,kk) = rmse(wk,w);

            if(corrTable(ii,jj,kk) > maxCorr)
                maxCorr = corrTable(ii,jj,kk);
                w_maxR = wk;
                %w_maxR = w_maxR .* (w_maxR > std(w_maxR)/2); % May remove
                %later
                ii1 = ii;
                jj1 = jj;
                kk1 = kk;
            end
            % if( g_dist(ii,jj,kk) < ming_dist )
            %     ming_dist = g_dist(ii,jj,kk);
            %     w_min_g_dist = wk;
            %     ll1 = ii;
            %     mm1 = jj;
            %     nn1 = kk;
            % end
            if( kendall_array(ii,jj,kk) > kendall_val )
                kendall_val = kendall_array(ii,jj,kk);
                w_kendall = wk;
                ll1 = ii;
                mm1 = jj;
                nn1 = kk;
            end
            if( mse_array(ii,jj,kk) < mse_val)
                mse_val = mse_array(ii,jj,kk);
                w_mse = wk;
                oo1 = ii;
                pp1 = jj;
                qq1 = kk;
            end
        end
    end
end

maxCorrArray = [maxCorr a(ii1) b(jj1) c(kk1)];
%ming_dist = [ming_dist a(ll1) b(mm1) c(nn1)];
kendall_valArray = [kendall_val a(ll1) b(mm1) c(nn1)];
mse_valArray = [mse_val a(oo1) b(pp1) c(qq1)];

w2d = squareform(w_maxR);
%w_min_g_dist2 = squareform(w_min_g_dist);
w_kendall = squareform(w_kendall);
w_mse = squareform(w_mse);

disp(subj)

%w_maxR = w_maxR .* (w_maxR > std(w_maxR)/5);

%[a1,b1,c1] = meshgrid(a,b,c);
%xslice = a(ii1);
%yslice = b(jj1);
%zslice = c(kk1);
%figure; slice(a1,b1,c1,corrTable,xslice,yslice,zslice);
%xlabel('\alpha' , 'FontSize' , 18); ylabel('\beta' , 'FontSize' , 18); zlabel('\gamma' , 'FontSize' , 18);
%colorbar

% if flagsin.fig == 1
%     figure;
%     sgtitle('SC recovery from time series, L1/L2' , 'FontSize' , 18);
% 
%     subplot(311);
%     plot(w);
%     title(['Vectorized SC, subject ' num2str(subj)]);
%     subplot(312);
%     plot(w_maxR);
%     title(['Est. SC via Pearson R = ' num2str(maxCorr, '%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ii1)) ', ' num2str(b(jj1)) ', ' num2str(c(kk1))]);
%     subplot(313);
%     plot(w_min_g_dist);
%     title(['Estimated SC, g_{dist} = ' num2str(ming_dist,'%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ll1)) ', ' num2str(b(mm1)) ', ' num2str(c(nn1))]);
% end
%corrBinary = corr(w_0,(w_maxR > 0));

if flagsin.fig == 1
    figure;
    %sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
    %subplot(2,2,1);
    imagesc( SCMtx / norm(SCMtx , 'fro'));
    %title(['True SC, subj ' num2str(subj)] )
    title( ['Empirical SC, subject ' num2str(subj)] );
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_' group '_subj_' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_ctrl_subj' num2str(subj)] , 'tif');
    saveas(gcf , [saveFigs filesep 'SC_' group '_subj_' num2str(subj)] , 'fig');
    clf

    figure;
    %subplot(2,2,2);
    imagesc( w2d / norm(w2d , 'fro') );
    title(['Est. SC, R = ' num2str(maxCorr,'%.2f') ', subj ' num2str(subj)]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj_' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_sub_j' num2str(subj)] , 'fig');
    clf

    % figure;
    % %subplot(2,2,3)
    % imagesc( w_min_g_dist2 / norm(w_min_g_dist2 , 'fro'));
    % title(['Riemann dist est. SC ' num2str(ming_dist,'%.2f') ', subj ' num2str(subj) ]);
    % axis square;
    % axis off;
    % saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj' num2str(subj)] , 'jpg');
    % %saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    % saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj' num2str(subj)] , 'fig');
    % clf

    figure;
    %subplot(2,2,3)
    imagesc( w_kendall / norm(w_kendall , 'fro'));
    title(['Kendall tau est. SC ' num2str(ming_dist,'%.2f') ', subj ' num2str(subj) ]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_kendall_' group '_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_kendall_' group '_subj' num2str(subj)] , 'fig');
    clf

    figure;
    %subplot(2,2,4)
    imagesc( w_mse / norm(w_mse , 'fro'));
    title(['MSE est. SC ' num2str(mve_val,'%.2f') ', subj ' num2str(subj)]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_mse_' group '_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_mse_' group '_subj' num2str(subj)] , 'tif');
    saveas(gcf , [saveFigs filesep 'SC_est_ts_FC_L1L2_mse_' group '_subj' num2str(subj)] , 'fig');
    clf
end

% figure;
% sgtitle('SC estimated via L1/L2, Riemann dist' , 'FontSize' , 15)
% subplot(2,2,1)
% imagesc( w_min_g_dist )
% title(['Riemann dist est. SC ' num2str(ming_dist)])
% axis square
% 
% subplot(2,2,2)
% imagesc(SCMtx);
% title(['True SC, subj ' num2str(subj)] )
% axis square
% 
% subplot(2,2,3)
% imagesc((w_min_g_dist > 0) .* (SCMtx > 0))
% title('Intersection')
% axis square
% 
% subplot(2,2,4)
% imagesc(((SCMtx > 0) .* (w_min_g_dist > 0)) - (w_min_g_dist > 0))
% title('Elements of est. SC are subset of true SC')
% axis square

% figure;
% sgtitle('SC estimated via L1/L2, R' , 'FontSize' , 15)
% subplot(2,2,1)
% imagesc( w2d )
% title(['Est. SC, R = ' num2str(maxCorr)])
% axis square
% 
% subplot(2,2,2)
% imagesc( SCMtx );
% title(['True SC, subj ' num2str(subj)] )
% axis square
% 
% subplot(2,2,3)
% imagesc((w2d > 0) .* (SCMtx > 0))
% title('Intersection')
% axis square
% 
% subplot(2,2,4)
% imagesc(( (SCMtx > 0) .* (w2d > 0) ) - (w2d > 0))
% title('Elements of est. SC are subset of true SC')
% axis square

%%

% if flagsin.fig == 1
%     figure;
%     %sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
%     subplot(1,3,1);
%     imagesc( SCMtx / norm(SCMtx));
%     %title(['True SC, subj ' num2str(subj)] )
%     title(['True SC'] )
%     axis square;
%     axis off
% 
%     subplot(1,3,2);
%     imagesc( w2d / norm(w2d) )
%     title(['Est. SC, R = ' num2str(maxCorr,'%.2f')])
%     axis square;
%     axis off
% 
%     subplot(1,3,3)
%     imagesc( w_min_g_dist2 / norm(w_min_g_dist2))
%     title(['Riemann dist est. SC ' num2str(ming_dist,'%.2f')])
%     axis square;
%     axis off
%     saveas(gcf , [saveFigs filesep 'SC_est_L1L2_subj' num2str(subj)] , 'jpg');
% end

%%

% if flagsin.fig == 1
%     figure;
%     sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
%     subplot(2,3,1);
%     imagesc( SCMtx / norm(SCMtx));
%     title(['True SC, subj ' num2str(subj)] )
%     axis square;
% 
%     subplot(2,3,2);
%     imagesc( w2d / norm(w2d) )
%     title(['Est. SC, R = ' num2str(maxCorr,'%.2f')])
%     axis square;
% 
%     subplot(2,3,3)
%     imagesc( w_min_g_dist2 / norm(w_min_g_dist2))
%     title(['Riemann dist est. SC ' num2str(ming_dist,'%.2f')])
%     axis square;
% 
%     subplot(2,3,4);
%     imagesc(squareform(w_0));
%     %colorbar;
%     title(['w_0 initial matrix, noise ' num2str(wstdn,'%.2f')]);
%     axis square;
% 
%     subplot(2,3,5)
%     imagesc( logical(w2d > 0) .* logical(SCMtx > 0) )
%     title('R intersection')
%     axis square
% 
%     subplot(2,3,6)
%     imagesc( logical(w_min_g_dist2 > 0) .* logical(SCMtx > 0) )
%     title('Riemann intersection')
%     axis square
% end
