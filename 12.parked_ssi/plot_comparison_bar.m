%% dash line plot
figurewidth = 18; %cm
f = figure('outerPosition',[10 10 figurewidth figurewidth*1]*29.1);
% subplot(2,2,1);subplot(2,2,2);
for i = 3:4
    
    if i==3
        subplot(2,2,i)
        [M1,I1] = max(pdfm784);
        [M2,I2] = max(arm784);
        [M3,I3] = max(acfm784);
        mn = [M1;M2;M3]';
        sd = [pdfs784(I1(1),1),pdfs784(I1(2),2),pdfs784(I1(3),3),pdfs784(I1(4),4);
            ars784(I2(1),1),ars784(I2(2),2),ars784(I2(3),3),ars784(I2(4),4);
            acfs784(I3(1),1),acfs784(I3(2),2),acfs784(I3(3),3),acfs784(I3(4),4)]';
    elseif i==4
        subplot(2,2,i)
        [M1,I1] = max(pdfm768);
        [M2,I2] = max(arm768);
        [M3,I3] = max(acfm768);
        mn = [M1;M2;M3]';
        sd = [pdfs768(I1(1),1),pdfs768(I1(2),2),pdfs768(I1(3),3),pdfs768(I1(4),4);
            ars768(I2(1),1),ars768(I2(2),2),ars768(I2(3),3),ars768(I2(4),4);
            acfs768(I3(1),1),acfs768(I3(2),2),acfs768(I3(3),3),acfs768(I3(4),4)]';
    end
    model_series = mn;
    model_error = sd;
    plot(1:4,model_series(:,1),'--o',1:4,model_series(:,2),'--o',...
        1:4,model_series(:,3),...
        '--o','MarkerSize',2,'linewidth',1);hold on
    nbars = size(model_series, 2);
    for j = 1:nbars
        if j == 1
            errorbar(1:4, model_series(:,j), model_error(:,j),...
                'color',[0, 0.4470, 0.7410],'linestyle', 'none','linewidth',1);
        elseif j == 2
            errorbar(1:4, model_series(:,j), model_error(:,j),...
                'color',[0.8500, 0.3250, 0.0980]	,'linestyle', 'none','linewidth',1);
        elseif j == 3
            errorbar(1:4, model_series(:,j), model_error(:,j),...
                'color',[0.9290, 0.6940, 0.1250],'linestyle', 'none','linewidth',1);
        end
    end
%     if i==3
%         title('(a) Damage level: 2\% stiffness reduction');
%     elseif i==4
%         title('(b) Damage level: 2\% stiffness reduction');
%     end
    ylim([0.5,1.02])
    yticks(0.5:0.1:1)
%     ,'Orientation','horizontal'
    xlim([0.9,4.1])
    xticks(1:4)
    xticklabels({'0%','15%','30%','45%'})
    xlabel('Noise level')
    ylabel('ROC-AUC')
%     if 
%     sgtitle(['(a) Damage level 1: $2\%$ stiffness reduction',...
% '(b) Damage level 2: $4\%$ stiffness reduction'],...
% 'interpreter','latex','location','south')
end

% % add a bit space to the figure
% fig = gcf;
% fig.Position(4) = fig.Position(4) + 250;
% % add legend
% Lgnd = legend('show');
% Lgnd.Position(1) = 0.01;
% Lgnd.Position(2) = 0.4;
% axes off
legend1 = legend('PDF-based','AR-based','ACF-based',...
        'location','northoutside','Orientation','horizontal');
set(legend1,...
    'Position',[0.203457439622135 0.486760261459986 0.640514191583539 0.0384382081997778],...
    'Orientation','horizontal',...
    'FontSize',7);
% sgtitle(['(a) Damage level 1: $2\%$ stiffness reduction',...
% '(b) Damage level 2: $4\%$ stiffness reduction'],...
% 'interpreter','latex','location','south')
set(findall(gcf,'-property','FontSize'),'FontSize',7)


exportgraphics(f,'fig7.comparison_errorbar.eps','Resolution',1000)

