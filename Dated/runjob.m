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

subj=1;
copyDir = '/Users/farrasabdelnour/Dropbox/Research//Projects'; 
FCdir = dir([copyDir filesep 'NYUEp_SC_to_FC' filesep 'results' filesep group corrType 'C' filesep '*_FC.mat']);
SCdir = dir([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep FCdir(subj).name(10:14) '*sc.csv']);


%dataLoc = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/Data/MEG/meg_individual_tseries_reordered.nc';
%finfo = ncinfo(dataLoc);
%vardata = ncread( dataLoc , values );

tsFCDir = dir('/Users/farrasabdelnour/Documents/Research/WCMC/Projects/NYUEpilepsyMtx/FunImgARWSDFC_AALTC_Ctrl/*AALTC.mat');
%SCMtx = csvread([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name],1,0);
SCMtx = readmatrix([copyDir filesep 'NYUEpilepsyMtx' filesep 'SC_Matrices' filesep SCdir(1).name]);

SCMtx = SCMtx(1:90,1:90);
stdSC = std(nonzeros(SCMtx(:)));
SCMtx = SCMtx .* (SCMtx > stdSC/5);

tsl = load([tsFCDir(subj).folder filesep tsFCDir(subj).name]);

%% Vectorize upper triangular
SCMtxt = SCMtx.';
m  = (1:size(SCMtxt,1)).' > (1:size(SCMtxt,2));
w  = SCMtxt(m);
clear m SCMtxt;


%% Vectorize upper triangle
%m  = (1:size(SCMtx,2)).' < (1:size(SCMtx,1));
%w  = SCMtx(m);
%clear m

w_0 = double(w > 0);

mtxFC = zeros(90,190);
for ii=1:90
    tt = eval(['tsl.AAL' num2str(ii,'%2.2d') 'TC']);
    mtxFC(ii,:) =  tt';
end

nSamples = size(mtxFC , 2);
X_noisy = mtxFC;

a = 0.001:0.05:2;
b = 0.001:0.05:2;
%c = 0.1:0.1:9;
maxCorr = 0;
corrTable = zeros(length(a) , length(b));
for ii = 1:length(a)
    for jj = 1:length(b)
        %for kk = 1:length(c)
        [wk, stat] = MM_gl(X_noisy, a(ii), b(jj),w_0,thresh,nSamples);
        corrTable(ii,jj) = corr(w,wk);
        if(abs(corrTable(ii,jj)) > maxCorr)
            maxCorr = corrTable(ii,jj);
            w_maxR = wk;
            ii1 = ii;
            jj1 = jj;
        end
        %end
    end
end

[maxR , Idx] = max( corrTable , [] , 'all');
[maxRx,maxRy] = find(corrTable == max( corrTable , [] , 'all'));
%maxRx = maxRx/length(a); maxRy = maxRy/length(b);

figure; subplot(211); plot(w); title(['Vectorized SC, subject ' num2str(subj)]); 
subplot(212); 
plot(w_maxR); title(['Estimated SC, R = ' num2str(maxR) ', \alpha, \beta = ' num2str(a(ii1)) ', ' num2str(b(jj1))]);

figure; 
imagesc( a , b, corrTable); cb = colorbar;
title(['Predicting SC from time series R, subj ' num2str(subj)] , 'FontSize' , 15)
xlabel('Parameter \alpha' , 'FontSize' , 15); ylabel('Parameter \beta' , 'FontSize' , 15)
axis square

thr = 0.999;
figure; 
imagesc(a,b,(corrTable > maxR*thr)); colorbar('Limits',cb.Limits);
title(['\alpha and \beta range with R \geq ' num2str(thr) 'Rmax, subj ' num2str(subj)] , 'FontSize' , 15);
xlabel('Parameter \alpha' , 'FontSize' , 15); ylabel('Parameter \beta' , 'FontSize' , 15)
axis square



