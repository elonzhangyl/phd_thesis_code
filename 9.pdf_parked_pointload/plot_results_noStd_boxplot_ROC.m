clear
load('result_pdf_KLDiv_k150_damage2and4.mat')
a2 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1:5);
a4 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,6:10);
load('result_pdf_KLDiv_k150_damage1and3.mat')
a1 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1:5);
a3 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,6:10);

b11 = cat(3,a1(:,:,1:25,:,:,:,1),a1(:,:,1:25,:,:,:,2),a1(:,:,1:25,:,:,:,3),a1(:,:,1:25,:,:,:,4),a1(:,:,1:25,:,:,:,5));
b21 = cat(3,a2(:,:,1:25,:,:,:,1),a2(:,:,1:25,:,:,:,2),a2(:,:,1:25,:,:,:,3),a2(:,:,1:25,:,:,:,4),a2(:,:,1:25,:,:,:,5));
b31 = cat(3,a3(:,:,1:25,:,:,:,1),a3(:,:,1:25,:,:,:,2),a3(:,:,1:25,:,:,:,3),a3(:,:,1:25,:,:,:,4),a3(:,:,1:25,:,:,:,5));
b41 = cat(3,a4(:,:,1:25,:,:,:,1),a4(:,:,1:25,:,:,:,2),a4(:,:,1:25,:,:,:,3),a4(:,:,1:25,:,:,:,4),a4(:,:,1:25,:,:,:,5));
b12 = cat(3,a1(:,:,26:50,:,:,:,1),a1(:,:,26:50,:,:,:,2),a1(:,:,26:50,:,:,:,3),a1(:,:,26:50,:,:,:,4),a1(:,:,26:50,:,:,:,5));
b22 = cat(3,a2(:,:,26:50,:,:,:,1),a2(:,:,26:50,:,:,:,2),a2(:,:,26:50,:,:,:,3),a2(:,:,26:50,:,:,:,4),a2(:,:,26:50,:,:,:,5));
b32 = cat(3,a3(:,:,26:50,:,:,:,1),a3(:,:,26:50,:,:,:,2),a3(:,:,26:50,:,:,:,3),a3(:,:,26:50,:,:,:,4),a3(:,:,26:50,:,:,:,5));
b42 = cat(3,a4(:,:,26:50,:,:,:,1),a4(:,:,26:50,:,:,:,2),a4(:,:,26:50,:,:,:,3),a4(:,:,26:50,:,:,:,4),a4(:,:,26:50,:,:,:,5));
b1 = cat(3,b11(:,:,1:100,:,:,:),b12(:,:,[1:44,46:101],:,:,:));
b2 = cat(3,b21(:,:,1:100,:,:,:),b22(:,:,1:100,:,:,:));
b3 = cat(3,b31(:,:,1:100,:,:,:),b32(:,:,1:100,:,:,:));
b4 = cat(3,b41(:,:,1:100,:,:,:),b42(:,:,1:100,:,:,:));
% b = cat(7,b1,b2,b3,b4);
b = cat(7,b2,b4);

%% boxplot-signle
for dd=1:2
figurewidth = 6; %cm
f = figure('Position',[10 10 figurewidth figurewidth*1]*36.36);
SCORE_avg = squeeze(mean(b,2));
for j=1:size(SCORE_avg,3)
    hold on
    scatter(j*ones(100,1),...
        squeeze(SCORE_avg(1,101:200,j,dd)),30,'x')
end
boxplot(squeeze(SCORE_avg(1,1:100,:,dd)),...
    'symbol','','Widths',0.3);hold on;
points = findobj(gcf, 'type', 'circle', 'Tag', 'Median');
set(points, 'Color', 'r');
xticks(1:4)
xticklabels({'0%','15%','30%','45%'})
ylim([0.8,3.5])
xlabel('Noise level')
ylabel('Damage index');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
tt = ['a','b','c','d'];
h1=findobj('Tag','Lower Whisker');
set(h1,'linestyle','-')
h2=findobj('Tag','Upper Whisker');
set(h2,'linestyle','-')
exportgraphics(f,['fig5',tt(dd),'.DI_Damage',num2str(dd),'.eps'],'Resolution',1000)
end

%% boxplot-all in one
figurewidth = 19; %cm
fig = figure('outerPosition',[1 1 figurewidth figurewidth*1.2]*36.36);
SCORE_avg = squeeze(mean(b,2));
t = tiledlayout(4,4,'TileSpacing', 'compact');
for ss=1:4
    for dd = 1:4
        nexttile
        for j=1:size(SCORE_avg,3)
            hold on
            scatter(j*ones(100,1),...
                squeeze(SCORE_avg(ss,101:200,j,dd)),10)
        end
        boxplot(squeeze(SCORE_avg(ss,1:100,:,dd)),...
            'symbol','','Widths',0.3);hold on;
        points = findobj(gcf, 'type', 'circle', 'Tag', 'Median');
        set(points, 'Color', 'r');
        xticks(1:4)
        xticklabels({'0%','15%','30%','45%'})
        set(gca,'fontsize',6)
        if ss==4
            ylim([0.8,15])
        else
            ylim([0.8,3.6])
        end
    end
end
xlabel(t, 'Noise level','FontSize',7)
ylabel(t, 'Damage index','FontSize',7)
% control dashed line whister
h1=findobj('Tag','Lower Whisker');
set(h1,'linestyle','-')
h2=findobj('Tag','Upper Whisker');
set(h2,'linestyle','-')
% set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(fig,['fig5','.DI_Damage','.eps'],'Resolution',1000)
% 
prob = zeros(4,4,4);
for ss=1:4
    for nl=1:4
        for dd=1:4
            pct = prctile(SCORE_avg(ss,1:100,nl,dd),95);
            prob(ss,nl,dd) = sum(SCORE_avg(ss,101:200,nl,dd)>pct);
        end
    end
end
pp = cat(1,prob(:,:,1)',prob(:,:,2)',prob(:,:,3)',prob(:,:,4)');

%% damage extent vs DI
% be = b(:,:,);

ss=1;
plot([squeeze(mean(SCORE_avg(ss,1:100,1,1),2));...
    squeeze(mean(SCORE_avg(ss,101:200,1,:),2))]);hold on
plot([squeeze(mean(SCORE_avg(ss,1:100,2,1),2));...
    squeeze(mean(SCORE_avg(ss,101:200,2,:),2))]);hold on
plot([squeeze(mean(SCORE_avg(ss,1:100,3,1),2));...
    squeeze(mean(SCORE_avg(ss,101:200,3,:),2))]);hold on
plot([squeeze(mean(SCORE_avg(ss,1:100,4,1),2));...
    squeeze(mean(SCORE_avg(ss,101:200,4,:),2))]);hold off
%% AUC - single
for dd=1:2 
NRT = 200;
tt = ['a','b','c','d'];
ss =1;
lstyle = {'-','--','-.',':'};
figurewidth = 7; %cm
fig = figure('outerPosition',[10 10 figurewidth figurewidth*1.2]*36.36);
SCORE_avg = squeeze(mean(b,2));
labels = (1:NRT)>(NRT/2);
% t = tiledlayout(4,4,'TileSpacing', 'compact');
    for nl = 1:size(SCORE_avg,3)% noise level
        [FPR,TPR,~,auc] = perfcurve(labels,SCORE_avg(ss,:,nl,dd),1);
        plot(FPR,TPR,'LineWidth',1,'linestyle',lstyle{nl});hold on
        if nl == 4
            hold off
        end
    end
xlabel('False positive rate','FontSize',7)
ylabel('True positive rate','FontSize',7)
legend(['  ','0% noise'],'15% noise','30% noise',...
    '45% noise','Orientation','vertical','location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(fig,['fig6',tt(dd),'.ROC_damage',num2str(dd),'.eps'],...
    'Resolution',1000)
end
%% AUC -all in one
NRT = 200;
lstyle = {'-','--','-.',':'};
figurewidth = 19; %cm
fig = figure('outerPosition',[1 1 figurewidth figurewidth*1]*29.1);
SCORE_avg = squeeze(mean(b,2));
labels = (1:NRT)>(NRT/2);
t = tiledlayout(4,4,'TileSpacing', 'compact');
% tlabel = nexttile([1,4]);
for ss=1:4
    for dd = 1:4
        nexttile        
        for nl = 1:size(SCORE_avg,3)% noise level
            [FPR,TPR,~,auc] = perfcurve(labels,SCORE_avg(ss,:,nl,dd),1);
            plot(FPR,TPR,'LineWidth',1,'linestyle',lstyle{nl});hold on
            if nl == 4
                hold off
            end
        end
    end
end
xlabel(t, 'False positive rate','FontSize',7)
ylabel(t, 'True positive rate','FontSize',7)
Lgnd = legend(nexttile(2),['  ','0% noise'],'15% noise','30% noise',...
    '45% noise','Orientation','horizontal','FontSize',7);
Lgnd.Layout.Tile = 'north';
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(fig,'fig6.ROC_damage.eps','Resolution',1000)



    