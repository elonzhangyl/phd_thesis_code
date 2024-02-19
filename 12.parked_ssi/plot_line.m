%% plot Dimensionality-ROC
figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.5]*29.1);
% axes('position',[0.1 0.1 0.9 0.9]);
xx = pca_order_list;
yy = Auc;
p = plot(xx,yy);
p.Color = [0, 0.4470, 0.7410];
p.LineStyle = '-';
p.LineWidth = 1;
p.Marker = 'o';
p.MarkerSize = 2.5;
p.MarkerEdgeColor = [0, 0.4470, 0.7410];
p.MarkerFaceColor = 'w';
% xlim([30,200]);
ylim([0.5 1]);
xticks(30:20:200);
xlabel('Dimensionality of the Hankel matrix, {\it{\fontname{Times} k}}');
ylabel('Average ROC-AUC');
% ytickformat('percentage')
% yline(0.05,'r');
% xlabel('Dimensionality of PDF, $k$','Interpreter', 'latex');
% ylabel('ROC-AUC');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig4.Dimensionality-ROC.eps','Resolution',1000)




