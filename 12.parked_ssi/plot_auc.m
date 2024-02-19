clear
dd = 8;
load('results_owt_KLIEP_PCA_0.24_pc14_sq130_de_dd10.mat')
a1 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
load('results_owt_KLIEP_PCA_0.36_pc14_sq130_de_dd10.mat')
a2 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
load('results_owt_KLIEP_PCA_0.48_pc14_sq130_de_dd10.mat')
a3 = SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd);
aa = cat(1,a1,a2,a3);
bb = squeeze(mean(aa,2));

NRT = 40;
labels = (1:NRT)>(NRT/2);
for dd=1:3
    [~,~,~,auc_temp] = perfcurve(labels,bb(dd,:),1);
    AUC(dd) = auc_temp; 
end
AUC

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.8]*29.1);
labels = (1:NRT)>(NRT/2);
dd_here = 1;
[FPR,TPR,~,auc] = perfcurve(labels,bb(1,:),1);
plot(FPR,TPR,'-.','linewidth',1);hold on
[FPR,TPR,~,auc] = perfcurve(labels,bb(2,:),1);
plot(FPR,TPR,'--','linewidth',1);hold on
[FPR,TPR,~,auc] = perfcurve(labels,bb(3,:),1);
plot(FPR,TPR,'-','linewidth',1);hold off
legend(['0.24 m damage'],...
    ['0.36 m damage'],...
    ['0.48 m damage'],...
    'Location','southeast');
xlabel('False positive rate');ylabel('True positive rate');
% title(['Sensor ',num2str(1)])
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig6a.ROC_owt.eps','Resolution',1000)
