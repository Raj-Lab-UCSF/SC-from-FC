


meanSC = mean(mtxSC,3);
stdSC = std(mtxSC,[], 3);

tt1 = meanSC(:);
pearsonCorr = zeros(1,26);
for jj = 1:size(mtxSC,3)
    tt = mtxSC(:,:,jj);
    tt = tt(:);
    
    pearsonCorr(jj) = corr(tt,tt1, 'type', 'Kendall');
end


for jj = 1:size(mtxSC,3)
    for ii = 1:size(mtxSC,3)

    tt1 = mtxSC(:,:,jj);
    tt2 = mtxSC(:,:,ii);    
    Pariwise_pearsonCorr(jj,ii) = corr(tt1(:),tt2(:), 'type', 'Kendall');
    end
end
figure; 
subplot(2,2,1); 
imagesc(meanSC); 
title('mean SC');

subplot(2,2,2); 
imagesc(stdSC);
title('std SC');

subplot(2,2,3); 
plot(pearsonCorr);
title('kendall tau')

subplot(2,2,4); 
imagesc(Pariwise_pearsonCorr);
title('Pairwise Kendall tau')

% % Now normalize each SC by sumsum
% 
% for jj = 1:size(mtxSC,3)
%     mtxSCn(:,:,jj) = mtxSC(:,:,jj)/sum(sum(mtxSC(:,:,jj),1),2);
% end
% meanSC = mean(mtxSCn,3);
% stdSC = std(mtxSCn,[], 3);
% 
% tt1 = meanSC(:);
% pearsonCorr = zeros(1,26);
% for jj = 1:size(mtxSCn,3)
%     tt = mtxSCn(:,:,jj);
%     tt = tt(:);
% 
%     pearsonCorr(jj) = corr(tt,tt1);
% end
% 
% 
% for jj = 1:size(mtxSCn,3)
%     for ii = 1:size(mtxSCn,3)
% 
%     tt1 = mtxSCn(:,:,jj);
%     tt2 = mtxSCn(:,:,ii);    
%     Pariwise_pearsonCorr(jj,ii) = corr(tt1(:),tt2(:));
%     end
% end
% figure; subplot(2,2,1); imagesc(meanSC); subplot(2,2,2); imagesc(stdSC)
% subplot(2,2,3); plot(pearsonCorr); 
% subplot(2,2,4); imagesc(Pariwise_pearsonCorr);
