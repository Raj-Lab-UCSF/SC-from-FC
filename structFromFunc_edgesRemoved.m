function [outw , meanSC, empiricalSC, x0_out , err] = structFromFunc_edgesRemoved( FC , mtxFC , meanSC , SC , p , subj, flagsin )

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

if ( nargin == 6)
    flagsin.fnctn = 'fmincon';
end

%numIter = 8;

vect = conv90altnodes2seq(90);
mtxFC = mtxFC(vect,:);
SC = conv90SPMalt2seqLR(SC);
FC = conv90SPMalt2seqLR(FC);
meanSC = conv90SPMalt2seqLR(meanSC);

mainFldr = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/';
saveFigs = [mainFldr 'fval_fconmin3'];

steps = logspace(-1,1,10);
fval = zeros(2,length(steps));

% Conservatively threshold meanSC
meanSC = meanSC .* (meanSC > std(meanSC,0,'all')/10);

%% Remove some inter hemispherical edges here:
%mask1 = zeros(size(meanSC));
%mask1(1:size(meanSC,1)/2 , size(meanSC,1)/2+1:end) = 1;
%mask1 = ~logical(mask1 + mask1'); 

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

%% Vectorize upper triangular mean and empirical SC

vecSC = meanSC;
vecSC = vectorizeMe(vecSC);
[~,mask] = sort(vecSC);
mask = (mask >= ceil(0.75 * length(mask) )); % was 0.85
%vecSC = vecSC(idx);
%vecSC = vectorizeMe(vecSC); % Inter-hemispheric now NaN

SC = vectorizeMe(SC); % empirical SC

%% Vectorize FC upper triangular 
d = vectorizeMe(FC);

%% Transpose time series
TS = mtxFC';

%% threshold x0, remove all negative elements of FC, threshold to std(x0)
%x0 = vectorizeMe(FC); %d; %w/2 + 0.25*std(w)*abs(randn(size(w,1),1)); % --> err = 0.2178, p=1, fmincon
x0 = FC; % * mean(meanSC(:))/mean(FC(:)); % x0 now has same mean as meanSC
x0 = vectorizeMe(x0);
x0 = x0 .* (abs(x0) > std(x0)); 
x0_out = squareform(x0);

%% Error, mean SC vs subject true SC
tempSC = vectorizeMe(meanSC);
err2 = norm( tempSC/norm(tempSC) - SC/norm(SC) );
clear tempSC;

% fMRI temporal resolution
params.TR = 2;

% In the following loop FC is first estimated from mean SC (fitSGM), 
% resulting estFC is vectorized (y), with L2 error between true (d) and y squared computed; 
% square root of w obtained and upper vectorized; and SC
% penalized as exp(-w).

kk = 1;
%nn = 1;

outerItr = 1;
err = zeros(outerItr,10); 

% lsqnonlin lower bound:
lb = -eps*ones(size(d));
ub = inf(size(x0)); 

tStart = tic;
[ theta_new , ~ , ~ , ~ , ~ ] = fitSGM( meanSC , TS , params ); % Here theta parameters are updated.
params.theta = theta_new;
tEnd = toc(tStart);

while kk <= outerItr

    %[ theta_new , ~ , ~ , ~ , ~ ] = fitSGM( updateSC , TS , params ); % Here theta parameters are updated.
    %params.theta = theta_new;

    ii = 1;
    % Test for jj = 2,3,4
    for jj = 1:10 %2:4 %1:length(steps) %0:numIter
        %for ll = steps % ll = 0:numIter

        lambda = [ii steps(jj)]; % ll];

        disp(['Lambda: ' num2str(steps(jj))])
        minErr = inf;

        if strcmp(flagsin.fnctn , 'lsqnonlin')
            options = optimoptions('lsqnonlin','Display','iter-detailed','MaxIterations', 24 ); %, 'MaxFunctionEvaluations', 5e4);
            [ w ,resnorm,residual,exitflag,output] = lsqnonlin( @(w)myfunc( y , w , d , vecSC  ,lambda , p) ,x0, lb, ub, options);
            fval(1,jj) = resnorm;

        elseif strcmp(flagsin.fnctn , 'fmincon')
            options = optimoptions('fmincon','Display','iter-detailed','MaxIterations', 24 , 'MaxFunctionEvaluations', 5e5 ,... 
                'PlotFcn',{@optimplotfval} , 'FunValCheck' , 'on' , 'HessianApproximation' , 'lbfgs');
            tfcon = tic;
            [ w , ~, ~,exitflag, ~] = fmincon( @(w)myfunc2( w , d , vecSC, params, TS, lambda, mask ), x0, [] , [], [], [], lb, ub, [], options);
            %[ w ,resnorm,residual,exitflag,output] = fmincon( @(w)myfunc2( y , w , d , vecSC, lambda, p ), x0, [] , [], [], [], lb, ub, [], options);
            fval(1,jj) = exitflag.bestfeasible.fval;
            tfconend = toc(tfcon)
        end


        % default algorithm: Algorithm','trust-region-reflective')
        % Other algorithms: 'interior-point', 'levenberg-marquardt'

        %% Comptue estimation error
        % SC: true SC, w: estimated subject SC
        err(kk,jj) = norm( w/norm(w) - SC/norm(SC) ); % L2 norm error, SC was vecSC

        if( err(kk,jj) <= minErr )
            minErr = err(kk,jj);
            outw = squareform(w); 
            minLambda = steps(jj);
        else
            outw = []; % No estimated SC error smaller than mean SC error
        end
        
        %fval(2,jj) = err(kk,jj);
        %kk = kk+1;
        % err needs to be smaller than err2sx
        disp(['Errors: ' num2str(minErr) num2str(err2)])
        %nn = nn+1;
        fig = gcf;
        figname = ['fmincon_itr' num2str(kk) '_x0_updated' num2str(subj) '_' num2str(jj)];
        saveas(gcf , [saveFigs filesep figname] , 'jpg');

    end

    %% Generate figure
    genfig(outw,squareform(SC),x0_out,meanSC,subj, minErr, minLambda, saveFigs);
    
    %% The following to be ignored for now!!
    %% Perform parameter udpate here
    %figure; imagesc(updateSC - squareform(w)); colorbar;
    %updateSC = squareform(w);
    %[ theta_new , estFC , ~ , ~ , ~ ] = fitSGM( updateSC , TS , params );

    % Update x0 here
    x0 = w/2; %/2 +  0.25*std(w)*abs(randn(size(w,1),1)); % x0 now not updated

    % Update the SC input to outer fitSGM...
    updateSC = squareform(w);

    %error1(nn) = tt;
    %nn = nn + 1;

    kk = kk+1;

    % Above to be ignored for now!
    
    disp(['kk = ' num2str(kk)])

end

%figure;
%plot(err);
%title(['Estimation err, subj ' num2str(subj) ', x0 updated'] , 'FontSize', 18);
%saveas(gcf , [saveFigs filesep 'error_x0_updated_subj_' num2str(subj)] , 'jpg');
%clf;

empiricalSC = squareform(SC);

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

function cost = myfunc2( w , d , vecSC, params, TS, lambda, mask )

% w: vectorized upper triangle mean SC to be estimated
% d: vectorized upper triangle FC
% vecSC: Upper triangle vectorized mean SC
% params: fitSGM parameters
% TS: fMRI time series
% lambda: fconmin weight parameter
% mask: mask selecting only the top 25% edges of mean SC

%% Threshold vectorized SC, w 
%stdw1 = std(w);
%w1 = squareform(w);
%w1 = w1 .* (w1 > stdw1/100);
%w1 = sum(logical(w1),2);
%w1 = sum(w1,2);

%updateSC = squareform(w); % Build estimated FC from latest SC estimate w
%tStart = tic;
[ ~ , estFC , ~ , ~ , ~ ] = fitSGM( squareform(w) , TS , params ); % here params.theta = [alpha,tau] are known
%tEnd = toc(tStart)
y = vectorizeMe(estFC);

%% Calculate difference between y and d
%rmsE = abs(y - d);

%% Calculate square root of the array w (equivalent to L1)
% p=1 corresponds to L2, p=2 corresponds to L1

% mask
%sqrtw = abs(w(mask) - vecSC(mask)); 

%% Exponential cost function, 
%expw = exp(-w1 / p );

%% x array:
temp = [lambda(1)* abs(y - d) ; lambda(2) * abs(w(mask) - vecSC(mask))]; 
%temp = [lambda(1) * abs(y - d)/norm(d) ; lambda(2) * abs(w(mask) - vecSC(mask))/norm(vecSC(mask))];

cost = norm(temp,2);

end


function x = myfunc3( y , w , d , meanSC, mask1, lambda, p )

% y: vectorized estimated FC from SGM
% w: vectorized upper diagonal mean SC to be estimated
% d: vectorized upper diagonal FC
% p: L1 (p=2) or L2 (p=1) imposed on the exp(w/p) term
% vecSC: Upper diagonal vectorized mean SC

%In that case you should add code in the myfun() that only uses non-NANs for cost function evaluation
% Use isnan() command
%

%% Threshold vectorized SC, w 
%stdw1 = std(w);
%w1 = squareform(w);
%w1 = w1 .* (w1 > stdw1/100);
%w1 = sum(logical(w1),2);
%w1 = sum(w1,2);

vecSC = meanSC .* mask1;
vecSC(~mask1) = NaN;
vecSC = vectorizeMe(vecSC); 

%% Calculate difference between y and d
rmsE = abs(y - d);

%% Calculate square root of the array w (equivalent to L1)
% p=1 corresponds to L2, p=2 corresponds to L1

vecSC = isnan(vecSC);
w = isnan(w);
sqrtw = abs(w - vecSC).^(1/p); 

%% Exponential cost function, 
%expw = exp(-w1 / p );


%% x Array:
x = [lambda(1)*rmsE ; lambda(2)*sqrtw]; % ; lambda(3)*expw];

x = norm(x,2);

end

function genfig(outw,empiricalSC,x0_out,meanSC,subj,minErr, lambda, saveFigs)

%
%
%
%

ftsize = 16;
figure;
    sgtitle(['Subj ' num2str(subj) ', lambda=' num2str(lambda) ', err='  num2str(minErr)] , 'FontSize' , 16);
    subplot(2,2,1);
    imagesc(outw); colorbar;
    title('Estimated SC' , 'FontSize',ftsize);

    subplot(2,2,2);
    imagesc(empiricalSC); colorbar;
    title('Empirical SC', 'FontSize',ftsize);

    subplot(2,2,3);
    imagesc(x0_out); colorbar;
    title('x0, initial guess' , 'FontSize',ftsize);

    subplot(2,2,4);
    imagesc(meanSC); colorbar;
    title('Mean SC', 'FontSize',ftsize);

    figname = ['SC_from_TS_'  num2str(subj)]; % '_' num2str(kk)];
    saveas(gcf , [saveFigs filesep figname] , 'jpg');

figure;
    sgtitle(['Subj ' num2str(subj) ', lambda=' num2str(lambda) ', err='  num2str(minErr) ', log'] , 'FontSize' , 16);
    subplot(2,2,1);
    imagesc(log(outw)); colorbar;
    title('Estimated SC, log' , 'FontSize',ftsize);

    subplot(2,2,2);
    imagesc(log(empiricalSC)); colorbar;
    title('Empirical SC, log', 'FontSize',ftsize);

    subplot(2,2,3);
    imagesc(x0_out); colorbar;
    title('x0, initial guess, log' , 'FontSize',ftsize);

    subplot(2,2,4);
    imagesc(log(meanSC)); colorbar;
    title('Mean SC, log', 'FontSize',ftsize);

    figname = ['SC_from_TS_log_'  num2str(subj)]; % '_' num2str(kk)];
    saveas(gcf , [saveFigs filesep figname] , 'jpg');


end