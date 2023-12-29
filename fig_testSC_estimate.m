function fig_testSC_estimate(outw,SC,meanSC)

figure; 
%up_clim = [-0.2 0.1];
subplot(3,3,1); 
imagesc( outw / norm( outw , 'fro') );  
axis square
set(gca,'XTick',[], 'YTick', [])
title('Estimated SC')
ylabel('Linear scale');

subplot(3,3,2); 
imagesc( SC / norm( SC , 'fro') ); 
axis square
set(gca,'XTick',[], 'YTick', [])
title('Empirical SC')

subplot(3,3,3); 
imagesc( meanSC / norm(meanSC, 'fro') );
axis square
set(gca,'XTick',[], 'YTick', [])
title('Mean SC')

%figure;
subplot(3,3,4);
%imagesc(log(outw .* (outw > std(outw,0,'all')/64) / norm(outw,'fro') + 1));
imagesc(log(outw / norm(outw,'fro')));
axis square
set(gca,'XTick',[], 'YTick', [])
%title('Estimated SC');
ylabel('Log scale');

subplot(3,3,5);
%imagesc(SC / norm(SC,'fro'));
%imagesc(log(SC .* (SC > std(SC,0,'all')/64) / norm(SC,'fro') + 1));
imagesc(log(SC / norm(SC,'fro') ));
axis square
set(gca,'XTick',[], 'YTick', [])
%title('Empirical SC');

subplot(3,3,6);
%imagesc(meanSC / norm(meanSC,'fro'));
%imagesc(log(meanSC .* (meanSC > std(meanSC,0,'all')/64) / norm(meanSC,'fro') + 1));
imagesc(log(meanSC / norm(meanSC,'fro') ));
set(gca,'XTick',[], 'YTick', [])
axis square
%title('Mean SC');

subplot(3,3,7);
imagesc(log(outw .* (outw > std(outw,0,'all')/64) / norm(outw,'fro')));
%imagesc(log(outw / norm(outw,'fro')));
axis square
set(gca,'XTick',[], 'YTick', [])

ylabel('Thresholded log')
subplot(3,3,8);
%imagesc(SC / norm(SC,'fro'));
imagesc(log(SC .* (SC > std(SC,0,'all')/64) / norm(SC,'fro')));
%imagesc(log(SC / norm(SC,'fro')));
axis square
set(gca,'XTick',[], 'YTick', [])
%title('Empirical SC');
subplot(3,3,9);
%imagesc(meanSC / norm(meanSC,'fro'));
imagesc(log(meanSC .* (meanSC > std(meanSC,0,'all')/64) / norm(meanSC,'fro')));
%imagesc(log(meanSC / norm(meanSC,'fro')));
set(gca,'XTick',[], 'YTick', [])
axis square
%title('Mean SC')





