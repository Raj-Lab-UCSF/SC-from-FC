function [maxCorr,ming_dist,Rcoor,gdistcoor] = runjob_l1_FC(subj,wstdn,flagsin)

% L1 pentalty with FC matrix replacing functional time series and *binarized
% FC* replacing (binarized) SC
%
% Use squareform to obtain matrix from vectorized w

group = 'Ctrl';
corrType = 'Pearson';
thresh = 0.001;
flagsin.SC = 1;
flagsin.FC = 0;
%
%subj=1;
%wstdn = 2^(-2);

%% Data locations
copyDir = '/Users/farrasabdelnour/Dropbox/Research/Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);
saveFigs = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/Figs';

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');
%SCMtx = csvread([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name],1,0);
%SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);

fsl = load([FCdir(subj).folder filesep FCdir(subj).name]);
FC = fsl.out.TrueFC;
SCMtx = fsl.out.StructConn;
stdSC = std(nonzeros(SCMtx(:)));
SCMtx = SCMtx .* (SCMtx > stdSC/5);

%tsl = load([tsFCDir(subj).folder filesep tsFCDir(subj).name]);

%% Vectorize upper triangular
%SCMtxt = SCMtx.';
%m  = (1:size(SCMtxt,1)).' > (1:size(SCMtxt,2));
%w  = SCMtxt(m);
%clear m SCMtxt fsl;

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


%% Vectorize upper triangle by stacking colmuns
% tf = tril(true(size(SCMtx,1)));
% SCMtx(tf) = NaN;
% idxnan = isnan(SCMtx);
% w = SCMtx(~idxnan);
% 
% clear tf

idxrnd = double(w_f == 0);
w_0 = double(w_f > 0) + (rand([size(w_f,1) 1])*wstdn) .* idxrnd;
%w_0 = zeros(length(w),1);

% Extract time series
%mtxFC = zeros(90,190);
%for ii=1:90
%    tt = eval(['tsl.AAL' num2str(ii,'%2.2d') 'TC']);
%    %tt = tsl.( 'AAL01TC' + num2str(ii,'%2.2d') + 'TC' );
%    mtxFC(ii,:) =  tt';
%end

nSamples = size(FC , 2);
%X_noisy = mtxFC;

a = (-3.001:0.05:3); % 
b = (-3.001:0.05:3); % 
%c = 0.001:0.02:1.5;
maxCorr = 0;
ming_dist = inf;
corrTable = zeros(length(a) , length(b)); % , length(c));
g_dist = zeros(length(a) , length(b));

for ii = 1:length(a)
    for jj = 1:length(b)
        %for kk = 1:length(c)
            [wk, ~] = MM_gl_l1_FC(FC, a(ii), b(jj),w_0,thresh,nSamples);
            wk = wk .* (wk > std(wk)/2); % May remove later
            corrTable(ii,jj) = corr(w,wk);
            g_dist(ii,jj) = riemann_dist( squareform(w) , squareform(wk) , flagsin );
            if(corrTable(ii,jj) > maxCorr)
                maxCorr = corrTable(ii,jj);
                w_maxR = wk;
                %alpha = a(ii);
                ii1 = ii;
                %beta = b(jj);
                jj1 = jj;
                Rcoor(1) = a(ii1);
                Rcoor(2) = b(jj1);
            end
            if( g_dist(ii,jj) < ming_dist && g_dist(ii,jj) > 0)
                ming_dist = g_dist(ii,jj);
                w_min_g_dist = wk;
                gdistcoor(1) = a(ii);
                gdistcoor(2) = b(jj);
            end
    end
end

%[maxR , Idx] = max( corrTable , [] , 'all');
%[maxRx,maxRy] = find(corrTable == max( corrTable , [] , 'all'));

if flagsin.fig == 1

figure; 
imagesc( [max(a) min(a)] , [max(b) min(b)], corrTable); colorbar;
title(['SC from FC, R = ' num2str(maxCorr) ', subj ' num2str(subj)] , 'FontSize' , 15)
xlabel('Parameter \alpha' , 'FontSize' , 15); ylabel('Parameter \beta' , 'FontSize' , 15)
axis square

figure; 
imagesc( [max(a) min(a)] , [max(b) min(b)], g_dist); 
oldcmap = colormap;
%colormap( flipud(oldcmap) );
colorbar;
title(['SC from FC, gdist = ' num2str(ming_dist) ', subj ' num2str(subj)] , 'FontSize' , 15)
xlabel('Parameter \alpha' , 'FontSize' , 15); ylabel('Parameter \beta' , 'FontSize' , 15)
axis square

% thr = 0.999;
% figure; 
% imagesc([max(a) min(a)], [max(b) min(b)], (corrTable > maxCorr*thr));
% title(['\alpha and \beta range with R \geq ' num2str(thr) 'Rmax, subj ' num2str(subj)] , 'FontSize' , 15);
% xlabel('Parameter \alpha' , 'FontSize' , 15); ylabel('Parameter \beta' , 'FontSize' , 15)
% axis square

%figure; subplot(211); plot(w); title(['Vectorized SC, subject ' num2str(subj)]); subplot(212); plot(w_maxR); 
%title(['Estimated SC, R = ' num2str(maxCorr) ', \alpha, \beta = ' num2str(a(ii1)) ', ' num2str(b(jj1))]);

%figure; subplot(211); plot(w); title(['Vectorized SC, subject ' num2str(subj)]); subplot(212); plot(w_min_g_dist); 
%title(['Estimated SC, g_dist = ' num2str(ming_dist) ', \alpha, \beta = ' num2str(a(kk1)) ', ' num2str(b(ll1))]);

w2d = squareform(w_maxR);
w_min_g_dist2 = squareform(w_min_g_dist);

% figure;
% subplot(2,2,1)
% imagesc((w2d > 0))
% title(['binarized est. SC, R= ' num2str(corrBinary)])
% axis square
% 
% subplot(2,2,2)
% imagesc((SCMtx > 0));
% title(['True SC binarized, subj ' num2str(subj)] )
% axis square
% 
% subplot(2,2,3)
% imagesc((w2d > 0) .* (SCMtx > 0))
% title('Intersection')
% axis square
% 
% subplot(2,2,4)
% imagesc(((SCMtx > 0) .* (w2d > 0)) - (w2d > 0))
% title('Elements of est. SC are subset of true SC')
% axis square


figure; 
%sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
subplot(1,3,1);
imagesc( SCMtx / norm(SCMtx));
title(['True SC, subj ' num2str(subj)] )
axis square;
axis off

subplot(1,3,2);
imagesc( w2d / norm(w2d) )
title(['Est. SC, R = ' num2str(maxCorr,'%.2f') ', L1'])
axis square;
axis off

subplot(1,3,3)
imagesc( w_min_g_dist2 / norm(w_min_g_dist2))
title(['Riemann dist est. SC ' num2str(ming_dist,'%.2f') ', L1'])
axis square;
axis off
saveas(gcf , [saveFigs filesep 'SC_est_L1_subj' num2str(subj)] , 'jpg');

end


