function [outw, empiricalSC, noisySC, x0_out, err, err2] = SC_fMRI_MICA_noise(subj , dataFldr , sparcity)

% This function opens MICA data mat file, extracts fMRI time series, SC,
% and FC for each subject and runs fitSGM for each subject.
%
%% Input:
% Subject subj of MICA is processed. If no subject is supplied, all subjects
% are processed.
%
%% output
% outw: 84x84x50 array of estimated SC for all subjects
% empiricalSC: 84x84x850 array of all subjects' empirical SC matrices
% x0_out: 84x84x50 array of subjects' fmincon initial guess
% err: estimation error between true SC and estimated SC
% err2: Error, mean SC vs subject true SC
%
% SC is normalized to sum/sum

%dataFldr = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';

%dataFldr = '../../Data/fMRI_MICA'; % Data location

numSubj = 50;
numNodes = 84;

S = load([dataFldr filesep 'MICA_aparc_struct.mat'] , "MICA" );
S = S.MICA;

flagsin.fnctn = 'fminunc';

%% Find mean SC

%tt = zeros(numNodes,numNodes,numSubj);
%for ii=1:numSubj
%   tt(:,:,ii) = S(ii).SC;
%end

meanSC = load("/wynton/home/rajlab/fabdelnour/data/fMRI_DK/mean_TFC.mat");

meanSC = meanSC.mean_TFC;
stdSC = meanSC.std_TFC;

%meanSC = mean(tt,3);
%stdSC = std(tt , 0, 3);

%zerosSC = (stdSC );

noise = 500 * abs(randn(numNodes));

% Sparcify noise
Index = randperm(numel(noise), ceil(numel(noise) * sparcity));
tt = noise( Index );
B = zeros( size(noise));
B(Index) = tt;
spNoise = B;

clear tt B Index ;

meanSC = ( meanSC + spNoise ) / mean( meanSC , "all");
stdSC = stdSC / mean( meanSC , "all" );

noisySC = meanSC;

numNodes = size(S(1).SC,1);

x0_out = zeros(numNodes,numNodes,size(S,2));
outw = x0_out;
empiricalSC = x0_out;

%err = Inf;
%err2 = 0;
% 
if ( nargin == 2 ) % Process all subjects
    err = zeros(numNodes,size(S,2));
    err2 = err;

    for subj=1:size(S,2)
        mtxFC = S(subj).fMRI_TS;
        %mtxFC = mtxFC';
        SC = S(subj).SC / sum( S(subj).SC , "all");
        FC = S(subj).FC;
        [outw(:,:,subj) , ~, empiricalSC(:,:,subj), x0_out(:,:,subj) , err(:,subj), err2(:,subj)] = structFromFunc_MICA( FC , mtxFC , meanSC , stdSC , SC  , subj, flagsin );
    end

elseif (nargin == 3) % Process only subject subj
    %err = zeros(numNodes,1);

    mtxFC = S(subj).fMRI_TS;
    %mtxFC = mtxFC';
    SC = S(subj).SC / sum( S(subj).SC , "all");
    FC = S(subj).FC;
    [outw , ~, empiricalSC, x0_out , err, err2 ] = structFromFunc_MICA( FC , mtxFC , meanSC , stdSC, SC  , subj, flagsin );

end




