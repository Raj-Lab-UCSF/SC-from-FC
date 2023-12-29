function [outw , meanSC, empiricalSC, x0 , fval] = structFromFunc_old( FC , mtxFC , meanSC , SC , p , flagsin )

%
% Functions called:
% lsqnonlin
% fmincon
% fitSGM
%
% Obtain NYU data using load_nyu_data.m
%
%% Input:
% mtxFC: times series -- mean time series can be obtaiend from mean_nyu_mtxFC.m 
% meanSC: mean empirical SC for error metric
% FC: Subject specific FC
% SC: Subject specific SC
% p = 1, L2 on w - meanSC term; p=2, L1 on w - meanSC term
% flagsin.fnctn: select 'lsqnonlin' or 'fmincon' (default)
%
%% Output:
% outw: the estimated SC
% metrics: estimation error and parameters lambda
% error1: estimation error
% fval: 2x10 array containing estimation error and lsqnonlin/fmincon iteration parameters
%
% To note: lambda weight inputs to lsqnonlin are squared
%
% EVALUATE both the error of estimated SC vs mean SC, and error of
% estimated SC vs true SC
%
% Figures to generate:
%
% outw, x0, meanSC, empirical SC, 
%
% Subjects 2, 5, 10, 20 

if ( nargin == 5)
    flagsin.fnctn = 'fmincon';
end

numIter = 8;

% Lp norm on exp(./p) penalty term . p=2 corresponds to L1, p=1 corresponds
% to L2

tt = inf;

%error1 = zeros((numIter+1)^3,1);
%metrics = zeros((numIter+1)^3 , 3); %4);


steps = logspace(-1,1,10);
fval = zeros(2,length(steps));

% Conservatively threshold meanSC
meanSC = meanSC .* (meanSC > std(meanSC,0,'all')/10);

%% Remove some inter hemispherical edges here:
% mask1 = zeros(size(meanSC));
% mask1(1:size(meanSC,1)/2 , size(meanSC,1)/2+1:end) = 1;
% rndMask = (abs(randn(size(meanSC))) > 1.2);
% rndMask = rndMask .* mask1;
% rndMask = rndMask + rndMask';
% rndMask = ~logical(rndMask);
% 
% clear mask1
% 
% meanSC = meanSC .* rndMask;
% 
% clear rndMask;

%empiricalSC = SC;

%% Vectorize upper triangular mean and empirical SC
vecSC = vectorizeMe(meanSC); % mean SC
SC = vectorizeMe(SC); % empirical SC

%% Vectorize upper triangular emperical FC
d = vectorizeMe(FC);

%% Binarize and add noise to upper triangular of mean SC
stdSC = std(meanSC(:));

w = vectorizeMe(meanSC);

%% Transpose time series
TS = mtxFC';
% Initialize
x0 = w/2 + 0.25*std(w)*abs(randn(size(w,1),1)); % --> err = 0.2178, p=1, fmincon

%% Error, mean SC vs subject true SC
err2 = norm( vecSC(:)/norm(vecSC) - SC(:)/norm(SC) );

% fMRI temporal resolution
params.TR = 2;

% In the following loop FC is first estimated from mean SC (fitSGM), 
% resulting estFC is vectorized (y), with L2 error between true (d) and y squared computed; 
% square root of w obtained and upper vectorized; and SC
% penalized as exp(-w).

%kk = 1;
nn = 1;

% lsqnonlin lower bound:
lb = -eps*ones(size(w));
ub = max(w(:))*ones(size(x0)) + 0.1;

[ theta_new , estFC , ~ , ~ , ~ ] = fitSGM( meanSC , TS , params ); % change meanSC updateSC

% Vectorize SGM estimated FC
y = vectorizeMe(estFC);
    
%for ii=0:numIter
ii = 1;
% Test for jj = 2,3,4
for jj = 3 %1:length(steps) %0:numIter
    %for ll = steps % ll = 0:numIter

    lambda = [ii steps(jj)]; % ll];
    disp([ii steps(jj)])
    tt = inf;

    %theta_new,

    if strcmp(flagsin.fnctn , 'lsqnonlin')
        options = optimoptions('lsqnonlin','Display','iter-detailed','MaxIterations', 24 ); %, 'MaxFunctionEvaluations', 5e4);
        [ w ,resnorm,residual,exitflag,output] = lsqnonlin( @(w)myfunc( y , w , d , vecSC  ,lambda , p) ,x0, lb, ub, options);
        fval(1,jj) = resnorm;

    elseif strcmp(flagsin.fnctn , 'fmincon')
        options = optimoptions('fmincon','Display','iter-detailed','MaxIterations', 24 , 'MaxFunctionEvaluations', 5e5 , 'PlotFcn',{@optimplotfval});
        [ w ,resnorm,residual,exitflag,output] = fmincon( @(w)myfunc2( y , w , d , vecSC, lambda, p ), x0, [] , [], [], [], lb, ub, [], options);
        fval(1,jj) = exitflag.bestfeasible.fval;
    end

    % default algorithm: Algorithm','trust-region-reflective')
    % Other algorithms: 'interior-point', 'levenberg-marquardt'

    %% Comptue estimation error
    % SC: true SC, w: estimated subject SC
    err = norm( w(:)/norm(w) - SC(:)/norm(SC) ); % L2 norm error, SC was vecSC

    %err_array_meanSC(kk) = err;
    %error1(nn) = err;
    if( err <= tt )
        tt = err;
        outw = squareform(w); % Computationally better to use meanSC
        % metrics(nn , :) = [err lambda];
    else
        outw = meanSC; % No estimated SC error smaller then mean SC error
        %metrics = [];
    end
    fval(2,jj) = tt;
    %kk = kk+1;
    % err needs to be smaller than err2
    disp([tt err2])
    %nn = nn+1;

    %% Build figure
    %figure;
    %imagesc(updateSC / norm(updateSC,2)); colorbar;
    %title(['itr ' num2str(kk-1) ', err1=' num2str(err,'%.2f') ', err2= ' num2str(err2,'%.2f') ', \lambda=' num2str(lambda(1),'%.1f') ', ' num2str(lambda(2),'%.1f') ', ' num2str(lambda(3),'%.1f')], 'fontsize' , 18);
    %axis square;
    %axis off;
    %pause;
    %clf;
end
%error1(nn) = tt;
nn = nn + 1;

x0 = squareform(x0);
empiricalSC = squareform(SC);
end
%end
%disp([ii jj ll])
% end

%end

function x = myfunc( y , w , d , vecSC, lambda , p )

% y: vectorized estimated FC from SGM
% w: vectorized upper diagonal mean SC to be estimated
% d: vectorized upper diagonal FC
% p: L1 (p=2) or L2 (p=1) imposed on the exp(w/p) term
% vecSC: Upper diagonal vectorized mean SC

%% Threshold vectorized SC, w 
% stdw1 = std(w);
% w1 = squareform(w);
% w1 = w1 .* (w1 > stdw1/100);
% %w1 = sum(logical(w1),2);
% w1 = sum(w1,2);

%% Calculate difference between y and d
rmsE = abs(y - d);

%% Calculate square root of the array w (equivalent to L1)
% p=1 corresponds to L2, p=2 corresponds to L1
%sqrtw = sqrt(w);
sqrtw = abs(w - vecSC).^(1/p);

%% Exponential cost function, 
%expw = exp(-w1 / p );

% Also try 1/w1 as an alternative to exp()
%expw = 1./sqrt((w1 + eps)).^p;
% How about
%expw = 1./(1 + w1.^(2/p));

%% x Array:
x = [lambda(1)*rmsE ; lambda(2)*sqrtw]; % ; lambda(3)*expw];

end

function vecSC = vectorizeMe(mtx)

% mtx is the matrix to be vectorized
% Only upper triangle elements are kept. This function assumes symmetric
% matrix.
%

vecSC = mtx';
m = (1:size(mtx,1)).' > (1:size(mtx,2));
vecSC = vecSC(m);

end

function x = myfunc2( y , w , d , vecSC, lambda, p )

% y: vectorized estimated FC from SGM
% w: vectorized upper diagonal mean SC to be estimated
% d: vectorized upper diagonal FC
% p: L1 (p=2) or L2 (p=1) imposed on the exp(w/p) term
% vecSC: Upper diagonal vectorized mean SC

%% Threshold vectorized SC, w 
%stdw1 = std(w);
%w1 = squareform(w);
%w1 = w1 .* (w1 > stdw1/100);
%w1 = sum(logical(w1),2);
%w1 = sum(w1,2);

%% Calculate difference between y and d
rmsE = abs(y - d);

%% Calculate square root of the array w (equivalent to L1)
% p=1 corresponds to L2, p=2 corresponds to L1

sqrtw = abs(w - vecSC).^(1/p); 

%% Exponential cost function, 
%expw = exp(-w1 / p );


%% x Array:
x = [lambda(1)*rmsE ; lambda(2)*sqrtw]; % ; lambda(3)*expw];

x = norm(x,2);

end
