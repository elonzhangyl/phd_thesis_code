clear
NRT = 40;

%% magnolia
% load('results_owt_KLIEP_PCA_0.48_pc13_sq130_de_dd10_alpha0.032.mat')
load('results_owt_KLD_0.24_pc8_sq160_de_dd10.mat')
labels = (1:NRT)>(NRT/2);
avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,:),2));



%% main
DD = 10;
metrics_values = zeros(DD,4);
for dd = 1:DD
    DI = avg(:,dd);
    DI_und = DI(1:20,:);
    DI_dam = DI(21:40,:);
    threshold  = mean(DI_und(:,1))+1.729*std(DI_und(:,1));

    TN = sum(DI_und<threshold);
    FP = sum(DI_und>=threshold);
    TP = sum(DI_dam>=threshold);
    FN = sum(DI_dam<threshold);
    accuracy = (TN+TP)./(TN+TP+FN+FP);
    recall =  TP./(TP+FN);
    f1score = 2*TP./(2*TP+FN+FP);
    MIOU = 1/2*(TN./(TN+FN+FP)+TP./(TP+FN+FP));

    metrics_values(dd,:) = [accuracy',recall',f1score',MIOU'];
end
metrics_values_mean = mean(metrics_values,1);
    % metrics_values_mean is the metrics for each damage level


