function [] = Fig_R_MM_parameters(params , corr )
%
%
% params is an array where first column contains Pearson R, Kendall tau, or MSE 
% remaining columns contain alpha, beta, and gamma MM parameters
% params = ts_SC.maxCorr_L2L1_ts_SC;
%


figure;
p = plot(params , 'o' , ...
    "LineWidth", 2, ...
    "MarkerSize" , 10);
legend([ corr "log, \alpha" "L2, \gamma" "L1, \beta"] , "AutoUpdate","on", 'FontSize', 16);
xlabel('Subjects' , 'FontSize', 16);
ylabel('MM Parameters' , 'FontSize', 16);

%set(gca,'xticklabel',{'Est SC', 'Mean SC','log','L2','L1'} , 'fontsize' , 6);


