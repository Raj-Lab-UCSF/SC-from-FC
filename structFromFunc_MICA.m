function [outw , meanSC, empiricalSC, x0_out , minErr, err2] = structFromFunc_MICA( FC , mtxFC , meanSC , stdSC, SC , subj, flagsin )

%
%% Functions called:
% 
% fitSGM
% squareform 
%
%% Input:
% mtxFC: times series, time x region.  
% meanSC: mean empirical SC for error metric
% FC: Subject specific FC
% SC: Subject specific SC
% flagsin.fnctn: select 'lsqnonlin' or 'fmincon' (default)
% subj: subject
%
%% Output:
% outw: the estimated SC correesponding to smallest err over lambda values
% meanSC: mean of subjects' SC
% empiricalSC: True SC
% x0_out: fmincon initial guess
% err: estimation error between true SC and estimated SC
% err2: Error, mean SC vs subject true SC
%
% To note: lambda weight inputs to lsqnonlin are squared
%
% Figures to generate:
%
% outw, x0, meanSC, empirical SC
%


if ( nargin == 6)
    flagsin.fnctn = 'fminunc';
end

% Lp norm on exp(./p) penalty term . p=2 corresponds to L1, p=1 corresponds
% to L2

%dataFldr = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';
%mainFldr = '/wynton/home/rajlab/fabdelnour/jobs/SC_from_fitSGM';

%dataFldr = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/Data/fMRI';
%mainFldr = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC';

mainFldr = pwd;
saveFigs = [mainFldr filesep 'Figs'];

steps = logspace(-1,1,10);
fval = zeros(2,length(steps));

% Conservatively threshold meanSC
%meanSC = meanSC; % .* (meanSC > std(meanSC,0,'all')/10);

%% Vectorize upper triangular mean, std, and empirical SC
vecSC = vectorizeMe(meanSC); % mean SC
vec_stdSC = vectorizeMe(stdSC);
Hstr = double(logical( vecSC )); % HessPattern fminunc 

% Keep the top 50% of the mean SC elements
pct = 0.5;
mask = (vecSC >= vecSC(ceil((1 - pct) * length(vecSC)) )); 

SC = vectorizeMe(SC); % vectorize empirical SC
d = vectorizeMe(FC); % vectorize FC

% Two lines below are meant to assign small values to the zeros as we
% discussed.
zerosSC = ( vec_stdSC < mean(vec_stdSC)/1000);
vec_stdSC(zerosSC) = mean(vec_stdSC)/1000;

%% Transpose time series
%TS = mtxFC';

%% threshold FC, zero out all negative elements of FC
x0 = d; 
%x0 = x0 * mean(vecSC)/mean(d);
x0 = x0 .* ( x0 > 0 );
%x0 = abs( x0 );
x0_out = squareform(x0);

%% Error, mean SC vs subject true SC
err2 = norm( vecSC/norm(vecSC) - SC/norm(SC) );

% fMRI temporal resolution for MICA data
params.TR = 0.6;

% In the following loop FC is first estimated from mean SC (fitSGM), 
% resulting estFC is vectorized (y), with L2 error between true (d) and y squared computed; 
% square root of w obtained and upper vectorized; and SC
% penalized as exp(-w).

%kk = 1;
%nn = 1;

outerItr = 1;
err = zeros(outerItr,1);

% lsqnonlin/fmincon bounds:
lb = -eps*ones(size(d));
ub = inf(size(x0)); 

tStart = tic;
[theta_new , ~ , ~ , ~ , ~ ] = fitSGM( meanSC , mtxFC , params ); % Here theta parameters are not updated.
params.theta = theta_new;
tEnd = toc(tStart);

%while kk <= outerItr

    %[ theta_new , ~ , ~ , ~ , ~ ] = fitSGM( meanSC , TS , params ); % Here theta parameters are updated.
    %params.theta = theta_new;

    %ii = 1;
    
    %for jj = 1:length(steps) % Evaluate estimated SC for all lambdas

    jj = 4; % The fourth value of steps array only
        %for ll = steps % ll = 0:numIter

        %lambda = steps(jj); % / mean(vec_stdSC(:)); 
        lambda = 0.01;
        disp(['lambda: ' num2str(lambda)])
        minErr = inf;

        if strcmp(flagsin.fnctn , 'lsqnonlin')
            options = optimoptions('lsqnonlin','Display','iter-detailed','MaxIterations', 32 ); %, 'MaxFunctionEvaluations', 5e4);
            [ w ,resnorm,residual,exitflag,output] = lsqnonlin( @(w)myfunc( y , w , d , vecSC  ,lambda , p) ,x0, lb, ub, options);
            fval(1,jj) = resnorm;

        elseif strcmp(flagsin.fnctn , 'fmincon')
            % MaxIterations set to 4 just for debugging purposes. Was 24
            options = optimoptions('fmincon','Display','iter-detailed','MaxIterations', 24 , 'MaxFunctionEvaluations', 5e5 ,... 
                'PlotFcn',{@optimplotfval} , 'FunValCheck' , 'on' , 'HessianApproximation' , 'lbfgs');
            [ w ,~,~,exitflag,~] = fmincon( @(w)myfunc2( w , d , vecSC, params, mtxFC , lambda, mask ), x0, [] , [], [], [], lb, ub, [], options);
            %[ w ,resnorm,residual,exitflag,output] = fmincon( @(w)myfunc2( y , w , d , vecSC, lambda, p ), x0, [] , [], [], [], lb, ub, [], options);
            fval(1,jj) = exitflag.bestfeasible.fval;
        
        elseif strcmp(flagsin.fnctn , 'fminunc')
             options = optimoptions('fminunc','Display','iter-detailed','MaxIterations', 19 , 'MaxFunctionEvaluations', 5e5 ,... 
                'PlotFcn',{@optimplotfval} , 'FunValCheck' , 'on' , 'HessianApproximation' , {"lbfgs",3}); % MaxIteration: 24
             %options = optimoptions('fminunc', 'Algorithm','trust-region', 'Display', 'iter-detailed','MaxIterations', 24 , 'MaxFunctionEvaluations', 5e5 ,... 
             %   'PlotFcn',{@optimplotfval} , 'FunValCheck' , 'on',  'SpecifyObjectiveGradient',true,'HessPattern',Hstr);
             %options = optimoptions('fminunc', 'Algorithm','trust-region',   'SpecifyObjectiveGradient',true,'HessPattern',Hstr);
           [w,fval,exitflag,output,grad,hessian] = fminunc(@(w)myfunc2( w , d , vecSC, vec_stdSC, params, mtxFC , lambda, mask ),x0,options);
        end


        % default algorithm: Algorithm','trust-region-reflective')
        % Other algorithms: 'interior-point', 'levenberg-marquardt'

        %% Compute estimation error between true SC and estimated SC, w
        % SC: true SC, w: estimated subject SC
        err(jj) = norm( w/norm(w) - SC/norm(SC) ); 

        if( err(jj) <= minErr )
            minErr = err(jj);
            outw = squareform(w); 
            minLambda = steps(jj);
        else
            outw = []; % No estimated SC error smaller then mean SC error
        end
        fval(2,jj) = minErr;
        %kk = kk+1;
        % err needs to be smaller than err2
        disp(['Errors: ' num2str(minErr) ' ' num2str(err2)])
        %nn = nn+1;

    %end
    
    %fig = gcf;
    %figname = ['fmincon_itr' '_' num2str(subj) '_' num2str(jj)];
    %saveas(gcf , [saveFigs filesep figname] , 'jpg');
    
    %% Generate figure
    %genfig(outw,squareform(SC),x0_out,meanSC,subj, minErr, lambda, saveFigs);
    %deltafig(outw,squareform(SC),meanSC,stdSC,subj,minErr, lambda, saveFigs);

    %% Perform parameter udpate here
    %figure; imagesc(updateSC - squareform(w)); colorbar;
    %updateSC = squareform(w);
    %[ theta_new , estFC , ~ , ~ , ~ ] = fitSGM( updateSC , TS , params ); 
    
    % Update x0 here
    %x0 = w/2; %/2 +  0.25*std(w)*abs(randn(size(w,1),1)); % x0 now not updated
    
    % Update the SC input to outer fitSGM...
    %updateSC = squareform(w);

    %kk = kk+1;
    %disp(['kk = ' num2str(kk)])
%end
 %figure; 
 %plot(err);
 %title(['Estimation err, subj ' num2str(subj)] , 'FontSize', 18);
 %saveas(gcf , [saveFigs filesep 'error_subj_' num2str(subj)] , 'jpg');
 %clf;

empiricalSC = squareform(SC);
end

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
x = [rmsE ; lambda*sqrtw]; % ; lambda(3)*expw];

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

function cost = myfunc2( w , d , vecSC, vec_stdSC, params, TS, lambda, mask )

% w: vectorized upper triangle SC to be estimated
% d: vectorized upper triangle FC
% vecSC: Upper triangle vectorized mean SC
% params: fitSGM parameters
% TS: fMRI time series
% lambda: fconmin weight parameter
% mask: mask selecting only the top 25% edges of mean SC

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

%% cost array:
%temp = [abs(y - d) ; lambda * abs(w(mask) - vecSC(mask))]; 
temp = [abs(y - d) ; lambda * abs((w - vecSC) ./ vec_stdSC) .* vecSC]; 
%temp = [abs(y - d)/norm(d) ; lambda * abs(w(mask) - vecSC(mask))/norm(vecSC(mask))];

cost = norm(temp,2);

end

function genfig(outw,empiricalSC,x0_out,meanSC,subj,minErr, lambda, saveFigs)

%
%
%
%

ftsize = 16;
figure;
    sgtitle(['Subj ' num2str(subj) ', lambda=' num2str(lambda) ', err='  num2str(minErr)] , 'FontSize' , ftsize );
    subplot(2,2,1);
    imagesc(outw); 
    c = colorbar;
    c.Limits = [0 max(empiricalSC(:))];
    title('Estimated SC' , 'FontSize',ftsize);

    subplot(2,2,2);
    imagesc(empiricalSC);
    c = colorbar;
    c.Limits = [0 max(empiricalSC(:))];
    title('Empirical SC', 'FontSize',ftsize);

    subplot(2,2,3);
    imagesc(x0_out); colorbar;
    title('x0, initial guess' , 'FontSize',ftsize);

    subplot(2,2,4);
    imagesc(meanSC); 
    c = colorbar;
    c.Limits = [0 max(empiricalSC(:))];
    title('Mean SC', 'FontSize',ftsize);

    figname = ['SC_from_TS_'  num2str(subj)]; % '_' num2str(kk)];
    saveas(gcf , [saveFigs filesep figname] , 'jpg');

figure;
    sgtitle(['Subj ' num2str(subj) ', lambda=' num2str(lambda) ', err='  num2str(minErr) ', log'] , 'FontSize' , ftsize);
    subplot(2,2,1);
    imagesc(log(abs( outw ))); 
    c = colorbar;
    c.Limits = [min(log(empiricalSC(:) + eps)) max(log(empiricalSC(:)))];
    title('Estimated SC, log' , 'FontSize',ftsize);

    subplot(2,2,2);
    imagesc(log(empiricalSC)); 
    c = colorbar;
    c.Limits = [min(log(empiricalSC(:) + eps)) max(log(empiricalSC(:)))];
    title('Empirical SC, log', 'FontSize',ftsize);

    subplot(2,2,3);
    imagesc(x0_out); colorbar;
    title('x0, initial guess' , 'FontSize',ftsize);

    subplot(2,2,4);
    imagesc(log(meanSC));
    c = colorbar;
    c.Limits = [min(log(empiricalSC(:) + eps)) max(log(empiricalSC(:)))];
    title('Mean SC, log', 'FontSize',ftsize);

    figname = ['SC_from_TS_log_'  num2str(subj)]; % '_' num2str(kk)];
    saveas(gcf , [saveFigs filesep figname] , 'jpg');


end


function deltafig(outw,empiricalSC,meanSC,stdSC,subj,minErr, lambda, saveFigs)

%
%
%
%

ftsize = 16;
figure;
    sgtitle(['Subj ' num2str(subj) ', lambda=' num2str(lambda) ', err='  num2str(minErr) ', \Delta'] , 'FontSize' , ftsize );
    subplot(2,2,1);
    imagesc( abs(meanSC - outw) ); colorbar;
    %c = colorbar;
    %c.Limits = [0 max(empiricalSC(:))];
    title('\Delta_{SC}' , 'FontSize',ftsize);

    subplot(2,2,2);
    imagesc( abs(meanSC - empiricalSC) ); colorbar;
    %c = colorbar;
    %c.Limits = [0 max(empiricalSC(:))];
    title('\Delta_{est SC}', 'FontSize',ftsize);

    subplot(2,2,3);
    imagesc( stdSC ); colorbar;
    title( 'SC group std', 'FontSize',ftsize );

    subplot(2,2,4);
    imagesc( meanSC ./ stdSC ); colorbar;
    title( 'Map of SC mean/std ratio' , 'FontSize',ftsize )

    figname = ['delta_SC_from_TS_'  num2str(subj)]; % '_' num2str(kk)];
    saveas(gcf , [saveFigs filesep figname] , 'jpg');

end


