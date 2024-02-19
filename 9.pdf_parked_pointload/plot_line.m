%% plot Dimensionality-ROC
figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.5]*29.1);
% axes('position',[0.1 0.1 0.9 0.9]);
xx = [1 5:5:300];
yy = squeeze(mean(AUC,5));
p = plot(xx,yy);
p.Color = [0, 0.4470, 0.7410];
p.LineStyle = '-';
p.LineWidth = 1;
p.Marker = 'o';
p.MarkerSize = 2.5;
p.MarkerEdgeColor = [0, 0.4470, 0.7410];
p.MarkerFaceColor = 'w';
xlim([1,250]);
% ylim([0 0.12]);
xticks([1,50:50:300]);
xlabel('Dimensionality of PDF, {\it{\fontname{Times} k}}');
ylabel('Average ROC-AUC');
% ytickformat('percentage')
% yline(0.05,'r');
% xlabel('Dimensionality of PDF, $k$','Interpreter', 'latex');
% ylabel('ROC-AUC');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig4.Dimensionality-ROC.eps','Resolution',1000)




