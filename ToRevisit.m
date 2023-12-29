

%% Decompose mean SC

[meanU , meanv] = eig(meanSC);
tt = zeros(90,90);
for ii = 1:size(mtxSCAll,3)
    tt = tt + (meanU' * mtxSCAll(:,:,ii) * meanU);
end

figure; imagesc(tt/size(mtxSCAll,3)); colorbar;

tt = tt/size(mtxSCAll,3);

d = eig(tt);

figure; 
sgtitle('Projection eigvenvalues vs diagonal' , 'fontsize' , 22);
subplot(2,3,1);
plot(sort(d) , 'LineWidth',2);
title('Mean proj eigenvals, d' , 'fontsize' , 14);

subplot(2,3,2);
plot(diag(tt) , 'LineWidth' , 2);
title('Diag of proj matrix, t' , 'fontsize' , 14);

subplot(2,3,3);
plot(sort(d) - diag(tt));
title('d - t' , 'FontSize', 14);

subplot(2,3,4);
imagesc(meanU' * mtxSCAll(:,:,7) * meanU );
title('NYU subj 7');

subplot(2,3,5);
imagesc(meanU' * mtxSCAll(:,:,5) * meanU);
title('NYU subj 5');

subplot(2,3,6);
imagesc(tt);
title('Mean of projection SC');


% evalsmean: eigenvalues of mean SC
% evals1: eig val of subj 2
% proj1 = diag(v'*mtxSCAll(:,:,2)*v);

figure; plot(evalsmean); hold on; plot(evals1); plot(proj1);


%% Filterbank wlet Matlab implementation.
% implement Ben's idea from my [1 -1] [1 1] Haar perspective

fb = dwtfilterbank('SignalLength',90,'Wavelet','db1','level',3);
[psi , t] = wavelets(fb);


figure; imagesc(psi);

x = conv(d2,psi(1,:));
figure;
plot(x(38:140));
hold on;
x = conv(d2,psi(2,:));
plot(x(38:140))
x = conv(d2,psi(3,:));
plot(x(38:140));
x = conv(d2,psi(4,:));
plot(x(38:140));
hold off;


[c,l] = wavedec(d10,3,'db1');
approx = appcoef(c,l,'db1');
[cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
figure;
fgtsize = 14;
sgtitle('SC eigenvalues Haar decomposition, subj 10' , 'FontSize' , 20);

subplot(4,1,1)
plot(approx , 'LineWidth', 2)
title('Approximation Coefficients' , 'FontSize', ftsize)
subplot(4,1,2)
plot(-cd3 , 'LineWidth', 2)
title('Bandpass level 3 detail coefficients' , 'FontSize', ftsize)
subplot(4,1,3 , 'LineWidth' , 2)
plot(-cd2 , 'LineWidth', 2)
title('Bandpass level 2 detail coefficients' , 'FontSize', ftsize)
subplot(4,1,4)
plot(-cd1 , 'LineWidth' , 2)
title('Highpass detail coefficients' , 'FontSize', ftsize)

%% Lift wavelet implementation
% sig: signal

lsc = liftingScheme('Wavelet','db1');
[ca,cd] = lwt(d10,'LiftingScheme',lsc,'Level',3);




