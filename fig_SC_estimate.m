function [] = fig_SC_estimate(outw , empiricalSC , meanSC , x0 , subj , flagsin )

%
% outw: estimated SC
% empSC: empirical SC
% meanSC: mean SC of a group (NYU for now)
% x0: fmincon initial guess
%
% Normalize matrices by Frobenius norm?

mainFldr = '/Users/farrasabdelnour/Library/CloudStorage/Box-Box/Research/UCSF/git/SC-from-FC/';
saveFigs = [mainFldr 'Figs/'];

subj = num2str(subj);

sgftsize = 18;
ftsize = 14;
yftsize = 14;

figure; 


t=tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'tight';

title(t,{'SC estimate', ['subject '  subj]} ,'FontSize' , sgftsize , 'Fontweight' , 'bold')
%up_clim = [-0.2 0.1];
%subplot(2,4,1); 
nexttile
imagesc( outw / norm( outw , 'fro') );  
axis square
set(gca,'XTick',[], 'YTick', [])
title('Estimated SC' , 'FontSize', ftsize );
ylabel('Linear scale' , 'FontSize', yftsize );

%subplot(2,4,2);
nexttile;
imagesc( empiricalSC / norm( empiricalSC , 'fro') ); 
axis square
set(gca,'XTick',[], 'YTick', [])
title('Empirical' , 'FontSize', ftsize );

%subplot(2,4,3);
nexttile;
imagesc( meanSC / norm(meanSC, 'fro') );
axis image
set(gca,'XTick',[], 'YTick', [])
title('Mean' , 'FontSize', ftsize );

%subplot(2,4,4);
nexttile;
imagesc( x0 / norm(x0,'fro'));
axis image
set(gca,'XTick',[], 'YTick', [])
title('Initial guess' , 'FontSize', ftsize );
%ylabel('Log scale');

%subplot(2,4,5);
nexttile;
imagesc(log(outw / norm(outw,'fro') ));
axis image
set(gca,'XTick',[], 'YTick', [])
%title('Estimated SC, log' , 'FontSize', ftsize );
ylabel('Log scale' , 'FontSize', yftsize );

%subplot(2,4,6);
nexttile;
imagesc(log(empiricalSC / norm(empiricalSC,'fro') ));
set(gca,'XTick',[], 'YTick', [])
axis image
%title('Empirical, log' , 'FontSize', ftsize );

%subplot(2,4,7);
nexttile;
imagesc(log( meanSC / norm(meanSC,'fro')));
axis image
set(gca,'XTick',[], 'YTick', [])
%title('Mean SC, log' , 'FontSize', ftsize );

%subplot(2,4,8);
nexttile;
imagesc(log( x0 / norm(x0,'fro')));
axis image
set(gca,'XTick',[], 'YTick', [])
%title('Initial guess, log' , 'FontSize', ftsize );

if flagsin.save == 1
    saveas(gcf , [saveFigs filesep 'SC_estimate_' subj] , 'jpg');
    saveas(gcf , [saveFigs filesep 'SC_estimate_' subj] , 'epsc');
end

