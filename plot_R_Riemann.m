function [] = plot_R_Riemann(t,s, flagsin)

%
% t is for Pearson R
% s is for Riemannian distance
%

saveFigs = [pwd filesep 'Figs'];
[Rs,d] = sort(t);
figure;

plot(Rs , 'LineWidth', 3); % Plot R

yyaxis left
ylabel('Pearson R', 'FontSize', 17)
yyaxis right
plot(s(d) , 'LineWidth', 3) % Riemannian
ylabel('Riemannian distance' , 'FontSize',17)
title('Pearson R & Riemannian metric' , 'FontSize', 17)
xlabel('Subjects')
set(gca,'FontSize',16);

if(flagsin.save == 1)
saveas(gcf , [saveFigs filesep 'subjects_R_Riemann_Ctrl'] , 'jpg');
saveas(gcf , [saveFigs filesep 'subjects_R_Riemann_Ctrl'] , 'epsc');
end