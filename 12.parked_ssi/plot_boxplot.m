clear
dd = 7;
load('results_owt_KLIEP_PCA_0.24_pc13_sq130_de_dd10_alpha0.032.mat')
a1 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
load('results_owt_KLIEP_PCA_0.36_pc13_sq130_de_dd10_alpha0.032.mat')
a2 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
load('results_owt_KLIEP_PCA_0.48_pc13_sq130_de_dd10_alpha0.032.mat')
a3 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
aa = cat(1,a1,a2,a3);
bb = squeeze(mean(aa,2));

NRT = 40;
labels = (1:NRT)>(NRT/2);
for dd=1:3
    [~,~,~,auc_temp] = perfcurve(labels,bb(dd,:),1);
    auc(dd) = auc_temp; 
end


figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.8]*29.1);
boxplot(squeeze([bb(1,1:20)',bb(1,1:20)',bb(1,1:20)']),...
    'symbol','','Widths',0.3);hold on;
points = findobj(gcf, 'type', 'circle', 'Tag', 'Median');
set(points, 'Color', 'r');
for j=1:size(bb,1)
    hold on
    scatter(j*ones(20,1),...
        squeeze(bb(j,21:40)),10)
end
xticks(1:4)
xticklabels({'0.24','0.36','0.48'})
% ylim([-inf,max(SCORE_rep_errorbar(:))])
ylim([3.308e-4,3.375e-4])
% yticks(2:6)
xlabel('Scour depth (m)')
%     if i == 1
ylabel('Damage index');
set(findall(gcf,'-property','FontSize'),'FontSize',7)

%% control dashed line whister
h1=findobj('Tag','Lower Whisker');
set(h1,'linestyle','-')
h2=findobj('Tag','Upper Whisker');
set(h2,'linestyle','-')

exportgraphics(f,'fig5a.DI_owt.eps','Resolution',1000)
