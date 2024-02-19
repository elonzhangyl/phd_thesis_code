clear
y = [0.638 0.751 0.823];
x = [0.599 0.725 0.819];

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.4]*36.36);
p1 = plot(1:3,x,'--o','MarkerSize',3,'linewidth',1);
hold on
p2 = plot(1:3,y,'-+','MarkerSize',4,'linewidth',1);
hold off

% ylim([0.6,1])
% yticks(0.7:0.1:1)
xlim([0.8,3.2])
xticks(1:3)
ylim([0.55 0.85])
xticklabels({'0.24','0.36','0.48'})
xlabel('Scour depth (m)');
ylabel('ROC-AUC')

% legend1 = legend('Gaussian-PDF-based','PCA-DRE-based',...
%         'location','southeast','Orientation','vertical');
    
legend([p2, p1], ["PCA-DRE-based", "Gaussian-PDF-based"],...
        'location','southeast','Orientation','vertical')

% set(legend1,...
%     'Position',[0.203457439622135 0.486760261459986 0.640514191583539 0.0384382081997778],...
%     'Orientation','horizontal',...
%     'FontSize',7);
legend boxoff
set(findall(gcf,'-property','FontSize'),'FontSize',7)

exportgraphics(f,'fig7.results_comparison_owt.eps','Resolution',1000)

