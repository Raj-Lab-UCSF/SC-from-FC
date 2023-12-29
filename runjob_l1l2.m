function [maxCorrArray,kendall_valArray,mse_valArray,w2d,w_kendall,w_mse] = runjob_l1l2(subj , subjSC , meanSC , flagsin)

% Here we perform L1/L2 penalty for estimating SC
%
% Majorization-minimization algorithm operates on *time series* and *L1/L2*
% on *binarized SC*
%
% w2d: estimated SC from Pearson R
% w_min_g_dist2: estimated SC from Riemann metric
%

group = 'Ctrl';
corrType = 'Pearson';
thresh = 0.001;
flagsin.SC = 1;
flagsin.FC = 0;

%wstdn = 2^(-3);

copyDir = '/Users/farrasabdelnour/Dropbox/Research/Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);

saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');
%SCMtx = csvread([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name],1,0);
SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);
SCMtx = SCMtx(1:90,1:90);

%stdSC = std(nonzeros(SCMtx(:)));
%SCMtx = SCMtx .* (SCMtx > stdSC/5);

tsl = load([tsFCDir(subj).folder filesep tsFCDir(subj).name]);

%% Vectorize upper triangular
meanSC = meanSC.';
m  = (1:size(meanSC,1)).' > (1:size(meanSC,2));
w  = meanSC(m);
clear m meanSC;

w_0 = w;

%% Vectorize empirical SC upper triangle
subjSC = subjSC.';
m  = (1:size(subjSC,1)).' > (1:size(subjSC,2));
w  = subjSC(m);
clear subjSC;

%% Vectorize upper triangle
%m  = (1:size(SCMtx,2)).' < (1:size(SCMtx,1));
%w  = SCMtx(m);
%clear m

% SC (here w) is now binarized 
%idxrnd = double(w == 0);
%w_0 = double(w > 0) + (rand([size(w,1) 1])*wstdn) .* idxrnd; % we're now adding some disruption to the zero elements of w_0

mtxFC = zeros(90,190); % fMRI Time series 
for ii=1:90
    tt = eval(['tsl.AAL' num2str(ii,'%2.2d') 'TC']);
    mtxFC(ii,:) =  tt';
end

nSamples = size(mtxFC , 2);

%f = -2.001;
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
            [wk, stat] = MM_gl_l2l1(mtxFC, a(ii), b(jj),c(kk),w_0,thresh,nSamples);
            wk = wk .* (wk > std(wk)/2); % May remove later
            corrTable(ii,jj,kk) = corr(w,wk);
            %kendall_array(ii,jj,kk) = riemann_dist( squareform(w) , squareform(wk) , flagsin );
            kendall_array(ii,jj,kk) = corr( w , wk , "type", "Kendall");
            mse_array(ii,jj,kk) = rmse(wk,w);

            if(corrTable(ii,jj,kk) > maxCorr)
                maxCorr = corrTable(ii,jj,kk);
                w_maxR = wk;
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
%     subplot(221);
%     plot(w);
%     title(['Vectorized SC, subject ' num2str(subj)]);
%     subplot(222);
%     plot(w_maxR);
%     title(['Est. SC via Pearson R = ' num2str(maxCorr, '%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ii1)) ', ' num2str(b(jj1)) ', ' num2str(c(kk1))]);
%     subplot(223);
%     plot(w_min_g_dist);
%     title(['Estimated SC, g_{dist} = ' num2str(ming_dist,'%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ll1)) ', ' num2str(b(mm1)) ', ' num2str(c(nn1))]);
% end

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

if flagsin.fig == 1
    figure;
    %sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
    %subplot(2,2,1);
    imagesc( SCMtx / norm(SCMtx , 'fro'));
    %title(['True SC, subj ' num2str(subj)] )
    title( ['Empirical SC, subject ' num2str(subj)] );
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_ctrl_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_ctrl_subj' num2str(subj)] , 'tif');
    %saveas(gcf , [saveFigs filesep 'SC_ctrl_subj' num2str(subj)] , 'fig');
    %clf

    figure;
    %subplot(2,2,2);
    imagesc( w2d / norm(w2d , 'fro') );
    title(['Est. SC, L1L2, R = ' num2str(maxCorr,'%.2f') ', subj ' num2str(subj)]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'fig');
    %clf

    figure;
    %subplot(2,2,3)
    imagesc( w_min_g_dist2 / norm(w_min_g_dist2 , 'fro'));
    title(['Riemann dist est. SC, L1L2 ' num2str(ming_dist,'%.2f') ', subj ' num2str(subj) ]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'fig');
    %clf

    figure;
    %subplot(2,2,4)
    imagesc( w_mse / norm(w_mse , 'fro'));
    title(['MSE est. SC ' num2str(mve_val,'%.2f') ', subj ' num2str(subj)]);
    axis square;
    axis off;
    saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'jpg');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'tif');
    %saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_R_' group '_subj' num2str(subj)] , 'fig');
    %clf
end




