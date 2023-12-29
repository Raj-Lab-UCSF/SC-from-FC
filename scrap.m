figure;
    %sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
    subplot(1,3,1);
    imagesc( mtxSC(:,:,1) / norm(mtxSC(:,:,1)));
    %title(['True SC, subj ' num2str(subj)] )
    title(['True SC'] )
    axis square;
    axis off

    subplot(1,3,2);
    imagesc( estSC_R / norm(estSC_R) )
    title(['Est. SC, R = ' num2str(maxCorrL2L1_ts_SC(1,1),'%.2f')])
    axis square;
    axis off

    subplot(1,3,3)
    imagesc( estSC_Riemann / norm(estSC_Riemann))
    title(['Riemann dist est. SC ' num2str(ming_distL2L1_ts_SC(1,1),'%.2f')])
    axis square;
    axis off


figure;
    sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
    subplot(2,3,1);
    imagesc( mtxSC(:,:,1) / norm(mtxSC(:,:,1) , 'fro'));
    title(['True SC, subj ' num2str(subj)] )
    axis square;

    subplot(2,3,2);
    imagesc( estSC_R / norm(estSC_R,'fro') )
    title(['Est. SC, R = ' num2str(maxCorrL2L1_ts_SC(1,1),'%.2f')])
    axis square;

    subplot(2,3,3)
    imagesc( estSC_Riemann / norm(estSC_Riemann , 'fro'))
    title(['Riemann dist est. SC ' num2str(ming_distL2L1_ts_SC(1,1),'%.2f')])
    axis square;

    subplot(2,3,4);
    imagesc(zeros(90,90));
    %colorbar;
    title(['w_0 initial matrix, noise ' num2str(90,'%.2f')]);
    axis square;

    subplot(2,3,5)
    imagesc( logical(estSC_R > 0.0001) .* logical(mtxSC(:,:,1) > 0.0001) )
    title('R intersection')
    axis square

    subplot(2,3,6)
    imagesc( logical(estSC_Riemann > 0.0001) .* logical(mtxSC(:,:,1) > 0.0001) )
    title('Riemann intersection')
    axis square


    for ii=1:26
        figure;
        imagesc(ts_SC.estSC_L2L1_R(:,:,ii)) ; colorbar;
        pause;
        %close;
    end


figure; 
up_clim = [0 0.4];
subplot(2,2,1); 
imagesc(meanSC / norm(meanSC,1) , up_clim); colorbar; 
axis square
title('Mean empirical SC')

subplot(2,2,2); 
imagesc(mean(estSC_L2L1_Kendall,3) / norm(mean(estSC_L2L1_Kendall,3), 1) , up_clim); colorbar;
axis square
title('Mean Kendall SC')

subplot(2,2,3); 
imagesc(mean(estSC_L2L1_R,3) / norm(mean(estSC_L2L1_R,3), 1) , up_clim); colorbar;
axis square
title('Mean Pearson SC')

subplot(2,2,4); 
imagesc(mean(estSC_L2L1_MSE,3) / norm(mean(estSC_L2L1_MSE,3), 1) , up_clim); colorbar;
axis square
title('Mean MSE SC')

for subj = 1:5

figure;
    tt = mtxSC(:,:,subj);
    %sgtitle('SC reconstruction using Pearson & Riemann metrics' , 'FontSize' , 15);
    subplot(1,3,1);
    imagesc( tt / norm(tt , 'fro'));
    %title(['True SC, subj ' num2str(subj)] )
    title( ['Empirical SC, subject ' num2str(subj)] );
    axis square;
    axis off;
    
    tt = ts_SC.estSC_L2L1_R(:,:,subj);
    subplot(1,3,2);
    imagesc( tt / norm(tt , 'fro') );
    title(['Est. SC, R = ' num2str(ts_SC.maxCorr_L2L1_ts_SC(subj,1),'%.2f') ', subj ' num2str(subj)]);
    axis square;
    axis off;
    
    tt = ts_SC.estSC_L2L1_MSE(:,:,subj);
    subplot(1,3,3)
    imagesc( tt / norm( tt , 'fro'));
    title(['MSE SC, L1L2, subj ' num2str(subj) ]);
    axis square;
    axis off;
    
 saveas(gcf , [saveFigs filesep 'SC_est_ts_meanSC_L1L2_subj' num2str(subj)] , 'jpg');

end

clear tt




if flagsin.fig == 1
figure;
sgtitle('3D visual of error from time series and FC, subj 5' , 'FontSize' , 18);
[a1,b1,c1] = meshgrid(a,b,c);
xslice = a(ii1);
yslice = b(jj1);
zslice = c(kk1);
figure; slice(a1,b1,c1,corrTable,xslice,yslice,zslice);
xlabel('\alpha' , 'FontSize' , 18); ylabel('\beta' , 'FontSize' , 18); zlabel('\gamma' , 'FontSize' , 18);
colorbar

if(flagsin.save == 1)
    %saveas(gcf , [saveFigs filesep 'L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'epsc');
    saveas(gcf , [saveFigs filesep '3D_L1L2_FC_from_ts_and_FC_subj' num2str(subj)] , 'jpg');
end
end



figure;
subplot(1,3,1);
imagesc(log(outw .* (outw > std(outw,0,'all')/64) / norm(outw,'fro')));
%imagesc(log(outw / norm(outw,'fro')));
axis square
set(gca,'XTick',[], 'YTick', [])
title('Estimated SC');

subplot(1,3,2);
%imagesc(SC / norm(SC,'fro'));
imagesc(log(SC .* (SC > std(SC,0,'all')/64) / norm(SC,'fro')));
%imagesc(log(SC / norm(SC,'fro')));
axis square
set(gca,'XTick',[], 'YTick', [])
title('Empirical SC');

subplot(1,3,3);
%imagesc(meanSC / norm(meanSC,'fro'));
imagesc(log(meanSC .* (meanSC > std(meanSC,0,'all')/64) / norm(meanSC,'fro')));
%imagesc(log(meanSC / norm(meanSC,'fro')));
set(gca,'XTick',[], 'YTick', [])
axis square
title('Mean SC')


w = w .* (w > std(w,0,'all')/128);

