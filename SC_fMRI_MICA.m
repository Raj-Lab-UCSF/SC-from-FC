function [outw, empiricalSC, x0_out, err, err2] = SC_fMRI_MICA(subj , dataFldr )

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

%dataFldr = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';

%dataFldr = '../../Data/fMRI_MICA'; % Data location

S = load([dataFldr filesep 'MICA_aparc_struct.mat'] , "MICA" );
S = S.MICA;

flagsin.fnctn = 'fminunc';

%% Find mean SC
%tt = zeros(size(S(1).SC));
%for ii = 1:size(S,2)
%    tt = tt + S(ii).SC;
%end
%meanSC = tt/size(S,2);
%clear tt;

tt = zeros(84,84,50);
for ii=1:50
    tt(:,:,ii) = S(ii).SC;
end

meanSC = mean(tt,3);
stdSC = std(tt , 0, 3);
%zerosSC = (stdSC );

clear tt;

numNodes = size(S(1).SC,1);

x0_out = zeros(numNodes,numNodes,size(S,2));
outw = x0_out;
empiricalSC = x0_out;

err = Inf;
err2 = 0;
% 
% if ( nargin == 1 ) % Process all subjects
%     err = zeros(numNodes,size(S,2));
%     err2 = err;
% 
%     for subj=1:size(S,2)
%         mtxFC = S(subj).fMRI_TS;
%         %mtxFC = mtxFC';
%         SC = S(subj).SC;
%         FC = S(subj).FC;
%         [outw(:,:,subj) , ~, empiricalSC(:,:,subj), x0_out(:,:,subj) , err(:,subj), err2(:,subj)] = structFromFunc_MICA( FC , mtxFC , meanSC , stdSC , SC  , subj, flagsin );
%     end
% 
% elseif (nargin == 2) % Process only subject subj
%     %err = zeros(numNodes,1);
% 
%     mtxFC = S(subj).fMRI_TS;
%     %mtxFC = mtxFC';
%     SC = S(subj).SC;
%     FC = S(subj).FC;
%     [outw , ~, empiricalSC, x0_out , err, err2 ] = structFromFunc_MICA( FC , mtxFC , meanSC , stdSC, SC  , subj, flagsin );
% 
% end




