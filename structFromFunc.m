function [outw , metrics, error1] = structFromFunc( FC , mtxFC , meanSC , SC)

%
% Functions called:
% lsqnonlin
% fitSGM
%
% Also consider optimoptions
%
% Input:
% mtxFC: times series -- mean time series can be obtaiend from mean_nyu_mtxFC.m 
% Obtain NYU data from function load_nyu_data.m
%
% meanSC: mean empirical SC for error metric
% FC: Subject specific FC
% SC: Subject specific SC
%
% To note: lambda weight inputs to lsqnonlin are squared
%
%
% EVALUATE both the error of estimated SC vs mean SC, and error of
% estimated SC vs true SC

%% Lambda lambda(1)*rmsE ; lambda(2)*sqrtw ; lambda(3)*expw
% lambda(1): rmsE
% lambda(2): L1 on w, estimated SC
% Lambda(3): L1 on w if p=2, L2 if p=1

numIter = 8;
% Lp norm on exp(./p) penalty term . p=2 corresponds to L1, p=1 corresponds
% to L2
p=1;
tt = inf;

error1 = zeros((numIter+1)^3,1);
metrics = zeros((numIter+1)^3 , 2); %4);

% steps = 0.1:0.1:10;
steps = logspace(-1,1,10);

% Conservatively threshold meanSC
meanSC = meanSC .* (meanSC > std(meanSC,0,'all')/10);

% Difference between SC and meanSC
D = SC - meanSC;
[Umean , ~] = eig(meanSC);

Beta = Umean' * D * Umean;
%Beta = vectorizeMe(Beta);

%% Remove some inter hemispherical edges here:
%mask1 = zeros(size(meanSC));
%mask1(1:size(meanSC,1)/2 , size(meanSC,1)/2+1:end) = 1;
%rndMask = (abs(randn(size(meanSC))) > 1.2);
%rndMask = rndMask .* mask1;
%rndMask = rndMask + rndMask';
%rndMask = ~logical(rndMask);
%clear mask1
%meanSC = meanSC .* rndMask;
%clear rndMask;

%% Vectorize upper trianngular mean and empirical SC
vecSC = vectorizeMe(meanSC); % mean SC
tSC = vectorizeMe(SC); % empirical SC

%% Vectorize upper triangular emperical FC
d = vectorizeMe(FC);

%w = vectorizeMe(meanSC);

%% Transpose time series
TS = mtxFC';
% Initialize
x0 = vecSC + 0.025*std(vecSC)*abs(randn(size(vecSC,1),1));

%% Error, mean SC vs subject true SC
err2 = norm( vecSC(:)/norm(vecSC , 'fro') - tSC(:)/norm(tSC , 'fro') );

% fMRI temporal resolution
params.TR = 2;

% In the following loop FC is first estimated from mean SC (fitSGM), 
% resulting estFC is vectorized (y), with L2 error between true (d) and y squared computed; 
% square root of w obtained and upper vectorized; and SC
% penalized as exp(-w).

%kk = 1;
nn = 1;

% lsqnonlin lower bound:
lb = -eps*ones(size(x0));
ub = max(x0(:))*ones(size(x0)) + 0.1;

[ theta_new , estFC , ~ , ~ , ~ ] = fitSGM( meanSC , TS , params ); 

% Vectorize SGM estimated FC
y = vectorizeMe(estFC);

%for ii=0:numIter
%ii = 1;
   for jj = steps %0:numIter
        %for ll = steps % ll = 0:numIter

            lambda = jj; % ll];
            disp(['\lambda_1 = ' jj])
            tt = inf;
            kk=1;
            while ( kk < 10 )
    
            %theta_new,
           
            % Initialize
            %x0 = w + 0.1*std(w)*abs(randn(size(w,1),1));
            
            % default algorithm: Algorithm','trust-region-reflective')
            % Algorithm: 'interior-point'
             %options = optimoptions('lsqnonlin','Display','iter','MaxIterations',24, 'MaxFunctionEvaluations', 5e4);
            % options.Algorithm = 'levenberg-marquardt'; % 
            options = optimoptions('fmincon','Display','iter','MaxIterations',24, 'MaxFunctionEvaluations', 5e4);

            %[ Beta ,resnorm,residual,exitflag,output] = lsqnonlin( @(Beta)myfunc( y , Beta , d  ,lambda , p) ,x0, lb, ub, options);
            [ Beta ,resnorm,residual,exitflag,output] = fmincon( @(Beta)myfunc2( y , Beta , d , lambda, p ), x0, [] , [], [], [], lb, ub, [], options);
    
   
            %% Comptue estimation error
            % SC: true SC, w: estimated subject SC
            estSC = Umean * squareform(Beta) * Umean' + meanSC;
            err = norm( estSC(:)/norm(estSC , 'fro') - SC(:)/norm(SC , 'fro') ); % L2 norm error, SC was vecSC

            %err_array_meanSC(kk) = err;
            %error1(nn) = err;
            if( err <= tt )
                tt = err; 
                outw = estSC; % Computationally better to use meanSC
                metrics(nn , :) = [err lambda];
            else
                outw = meanSC; % No estimated SC error smaller then mean SC error
                %metrics = [];
            end

            kk = kk+1;
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
            error1(nn) = tt;
            nn = nn + 1;
        end
    %end
    %disp([ii jj ll])
% end

end

function x = myfunc( y , d , Beta, lambda , p )

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

%% Calculate L1/L2 of Beta
%sqrtw = sqrt(w);
%sqrtw = abs(w - vecSC);

Beta = Beta.^(1/p);

%% Exponential cost function, 
%expw = exp(-w1 / p );

% Also try 1/w1 as an alternative to exp()
%expw = 1./sqrt((w1 + eps)).^p;
% How about
%expw = 1./(1 + w1.^(2/p));

%% x Array:
x = [ rmsE ; lambda * Beta]; % ; lambda(3)*expw];

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

function x = myfunc2( y  , d , Beta, lambda, p )

% y: vectorized estimated FC from SGM
% w: vectorized upper triangular mean SC to be estimated
% d: vectorized upper triangular FC
% p: L1 (p=2) or L2 (p=1) imposed on the exp(w/p) term
% vecSC: Upper diagonal vectorized mean SC


%% Calculate difference between y and d
rmsE = abs(y - d);

%% Calculate square root of the array Beta (equivalent to L1)
%sqrtw = sqrt(w);
%sqrtw = sqrt(abs(w - vecSC)); % L1 case
%sqrtw = abs(w - vecSC)^(1/p); % This is the L2 case

Beta = Beta.^(1/p);

%% x Array:
x = [rmsE ; lambda * Beta]; % ; lambda(3)*expw];

x = norm(x,2);

end
