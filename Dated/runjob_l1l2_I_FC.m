function [maxCorr,ming_dist] = runjob_l1l2_I_FC( subj , wstdn , flagsin )

%runjob

%
%
%
%
%
%
%

group = 'Ctrl';
corrType = 'Pearson';
thresh = 0.001;
flagsin.SC = 1;
flagsin.FC = 0;

%wstdn = 2^(-3);

%subj=1;
copyDir = '/Users/farrasabdelnour/Dropbox/Research/Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
%SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');
%SCMtx = csvread([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name],1,0);
%SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);

saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';

%SCMtx = SCMtx(1:88,1:88);
%stdSC = std(nonzeros(SCMtx(:)));
%SCMtx = SCMtx .* (SCMtx > stdSC/5);

% Load FC of subject subj
%tsl = load([tsFCDir(subj).folder filesep tsFCDir(subj).name]);
fsl = load([FCdir(subj).folder filesep FCdir(subj).name]);
FC = fsl.out.TrueFC;
SCMtx = fsl.out.StructConn;
stdSC = std(nonzeros(SCMtx(:)));
SCMtx = SCMtx .* (SCMtx > stdSC/5);

%% Vectorize SCM upper triangular
SCMtxt = SCMtx.';
m  = (1:size(SCMtxt,1)).' > (1:size(SCMtxt,2));
w  = SCMtxt(m);
clear m SCMtxt fsl;

%% Vectorize FC upper triangular
% FC = abs(FC);
% 
absFC = abs(FC);
stdFC = std(nonzeros(absFC(:)));
threshFC = absFC .* (absFC > 2*stdFC);
threshFC = threshFC.';
m  = (1:size(threshFC,1)).' > (1:size(threshFC,2));
w_f  = threshFC(m);
clear m;


%% Vectorize upper triangle
%m  = (1:size(SCMtx,2)).' < (1:size(SCMtx,1));
%w  = SCMtx(m);
%clear m

%w_0 = double(w > 0);
%idxrnd = double(w == 0);
%w_0 = double(w > 0) + (rand([size(w,1) 1])*wstdn) .* idxrnd; % we're now adding some disruption to the zero elements of w_0

% w_0 now is based on thresholded FC now
idxrnd = double(w_f == 0);
w_0 = double(w_f > 0) + (rand([size(w_f,1) 1])*wstdn) .* idxrnd; % we're now adding some disruption to the zero elements of w_0

%mtxFC = zeros(90,190);
%for ii=1:90
%    tt = eval(['tsl.AAL' num2str(ii,'%2.2d') 'TC']);
%    mtxFC(ii,:) =  tt';
%end

nSamples = size(FC , 2);
%X_noisy = mtxFC;

%f = 0.001;
f = 0.001;
a = (f:0.05:1);
b = (f:0.05:1);
c = (f:0.05:1);
maxCorr = 0;
ming_dist = inf;
corrTable = zeros(length(a) , length(b) , length(c)); 
g_dist = zeros(length(a) , length(b) , length(c));
for ii = 1:length(a)
    for jj = 1:length(b)
        for kk = 1:length(c)
            %[wk, stat] = MM_gl_l2l1(mtxFC, a(ii), b(jj),c(kk),w_0,thresh,nSamples);
            [wk, stat] = MM_gl_l2l1_I_FC(FC, a(ii), b(jj),c(kk),w_0,thresh,nSamples);
            wk = wk .* (wk > std(wk)/2); % May remove later
            corrTable(ii,jj,kk) = corr(w,wk);
            g_dist(ii,jj,kk) = riemann_dist( squareform(w) , squareform(wk) , flagsin );
            if(corrTable(ii,jj,kk) > maxCorr)
                maxCorr = corrTable(ii,jj,kk);
                w_maxR = wk;
                %w_maxR = w_maxR .* (w_maxR > std(w_maxR)/2); % May remove later
                alpha = a(ii);
                ii1 = ii;
                beta = b(jj);
                jj1 = jj;
                %kk1 = kk;
                gamma = c(kk);
                kk1 = kk;
            end
            if( g_dist(ii,jj,kk) < ming_dist )
                ming_dist = g_dist(ii,jj,kk);
                w_min_g_dist = wk;
                ll1 = ii;
                mm1 = jj;
                nn1 = kk;
            end
        end
    end
end

%w_maxR = w_maxR .* (w_maxR > std(w_maxR)/5);

if flagsin.fig == 1
figure;
sgtitle(['3D visual of L1/L2 from time series and FC, subj ' num2str(subj)] , 'FontSize' , 18);
[a1,b1,c1] = meshgrid(a,b,c);
xslice = a(ii1);
yslice = b(jj1);
zslice = c(kk1);
figure; slice(a1,b1,c1,corrTable,xslice,yslice,zslice);
xlabel('\alpha' , 'FontSize' , 18); ylabel('\beta' , 'FontSize' , 18); zlabel('\gamma' , 'FontSize' , 18);
colorbar

if(flagsin.save == 1)
    %saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'epsc');
    saveas(gcf , [saveFigs filesep '3D_L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'jpg');
end
end

%% Just a bunch of plot figures
% figure; 
% sgtitle(['L1/L2 SC from time series and FC, subj ' num2str(subj)] , 'FontSize' , 18);
% 
% subplot(3,1,1); 
% plot(w); 
% title(['Vectorized SC, subject ' num2str(subj)]); 
% subplot(3,1,2); 
% plot(w_maxR); 
% title(['Est. SC via Pearson R = ' num2str(maxCorr, '%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ii1)) ', ' num2str(b(jj1)) ', ' num2str(c(kk1))]);
% subplot(3,1,3); 
% plot(w_min_g_dist); 
% title(['Estimated SC, g_{dist} = ' num2str(ming_dist,'%.2f') ', \alpha, \beta, \gamma = ' num2str(a(ll1)) ', ' num2str(b(mm1)) ', ' num2str(c(nn1))]);
% 
% if(flagsin.save == 1)
%     %saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'epsc');
%     saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'jpg');
% end

%corrBinary = corr(w_0,(w_maxR > 0));

w2d = squareform(w_maxR);
w_min_g_dist2 = squareform(w_min_g_dist);

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

if flagsin.fig == 1
    figure;
    sgtitle(['L1/L2 SC from time series and FC, subj ' num2str(subj)] , 'FontSize' , 15);
    subplot(2,3,1);
    imagesc( SCMtx / norm(SCMtx));
    title(['True SC, subj ' num2str(subj)] )
    axis square;

    subplot(2,3,2);
    imagesc( w2d / norm(w2d) )
    title(['Est. SC, R = ' num2str(maxCorr,'%.2f')])
    axis square;

    subplot(2,3,3)
    imagesc( w_min_g_dist2 / norm(w_min_g_dist2))
    title(['Riemann dist est. SC ' num2str(ming_dist,'%.2f')])
    axis square;

    subplot(2,3,4);
    imagesc(squareform(w_0));
    %colorbar;
    title(['w_0 initial matrix, noise ' num2str(wstdn,'%.3f')]);
    axis square;

    subplot(2,3,5)
    imagesc( logical(w2d > 0) .* logical(SCMtx > 0) )
    title('R intersection')
    axis square

    subplot(2,3,6)
    imagesc( logical(w_min_g_dist2 > 0) .* logical(SCMtx > 0) )
    title('Riemann intersection')
    axis square

    if(flagsin.save == 1)
        %saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'epsc');
        saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'jpg');
    end
end


