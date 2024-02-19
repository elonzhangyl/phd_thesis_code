clear;clc;

tic;
%% Tuning parameter
sqenc_list = 205;
alpha_list = [0 0.15 0.30 0.45];
sensor_list = 1:4;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;
NRR = 25;
NRT = 50; 
DD = 10;
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
for dd = 1:DD
    dd
    file_ref = {'und_no_water_real_6000s_seed201To205.mat',...
        'und_no_water_real_6000s_seed211To215.mat',...
        'und_no_water_real_6000s_seed221To225.mat',...
        'und_no_water_real_6000s_seed231To235.mat',...
        'und_no_water_real_6000s_seed241To245.mat',...
        'und_no_water_real_6000s_seed201To205.mat',...
        'und_no_water_real_6000s_seed211To215.mat',...
        'und_no_water_real_6000s_seed221To225.mat',...
        'und_no_water_real_6000s_seed231To235.mat',...
        'und_no_water_real_6000s_seed241To245.mat'};
    load(string(file_ref(dd)));
    dataset_ref=dataset1;
    dataset_ref = reshape(dataset_ref,1200*20,25,4);
    
    file_test_1 = {'und_no_water_real_6000s_seed206To210.mat',...
        'und_no_water_real_6000s_seed216To220.mat',...
        'und_no_water_real_6000s_seed226To230.mat',...
        'und_no_water_real_6000s_seed236To240.mat',...
        'und_no_water_real_6000s_seed246To250.mat',...
        'und_no_water_real_6000s_seed206To210.mat',...
        'und_no_water_real_6000s_seed216To220.mat',...
        'und_no_water_real_6000s_seed226To230.mat',...
        'und_no_water_real_6000s_seed236To240.mat',...
        'und_no_water_real_6000s_seed246To250.mat'};
    load(string(file_test_1(dd)));
    dataset_test_1 = dataset1;
    dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);


    file_test_2 = {'dam_no_water_real_allmatrix_soil7.84_seed301To305.mat',...
        'dam_no_water_real_allmatrix_soil7.84_seed306To310.mat',...
        'dam_no_water_real_allmatrix_soil7.84_seed311To315.mat',...
        'dam_no_water_real_allmatrix_soil7.84_seed316To320.mat',...
        'dam_no_water_real_allmatrix_soil7.84_seed321To325.mat',...
        'dam_no_water_real_allmatrix_soil7.68_seed401To405.mat',...
        'dam_no_water_real_allmatrix_soil7.68_seed406To410.mat',...
        'dam_no_water_real_allmatrix_soil7.68_seed411To415.mat',...
        'dam_no_water_real_allmatrix_soil7.68_seed416To420.mat',...
        'dam_no_water_real_allmatrix_soil7.68_seed421To425.mat'};
    load(string(file_test_2(dd)));

    dataset_test_2 = dataset2;
    dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);
    dataset_test = cat(2,dataset_test_1,dataset_test_2);

  
    %% Damage Index
    SCORE_SS_NRR_NRT_SN_NL_SQ = zeros(SS,NRR,NRT,SN,NL,SQ);
    for sq = 1:SQ
        sqenc = sqenc_list(sq);
        SCORE_SS_NRR_NRT_SN_NL = zeros(SS,NRR,NRT,SN,NL);
        for nl = 1:NL % noise level
%             nl
            alpha =alpha_list(nl);
            SCORE_SS_NRR_NRT_SN = zeros(SS,NRR,NRT,SN);
            for sn = 1:SN % seed noise
                sn
                SCORE_SS_NRR_NRT = zeros(SS,NRR,NRT);
                parfor nrt=1:NRT % NUM RECORD TEST
                    nrt
                    Dataset2 = dataset_test(:,nrt,:); %choose record 
    %                 rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_SS_NRR = zeros(SS,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,nrr,:); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = zscore(Dataset10); %normlization
                        SCORE_SS = zeros(SS,1);
                        for ss = 1:SS
                            X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                            X2 = zeros(size(Dataset200,1)-sqenc,sqenc);
                            for j=1:(size(Dataset100,1)-sqenc)
                                X1(j,:) = Dataset100(j:j+sqenc-1,ss);
                            end
                            for j=1:(size(Dataset200,1)-sqenc)
                                X2(j,:) = Dataset200(j:j+sqenc-1,ss);
                            end
                        SCORE_SS(ss) = PEDiv(X1,X2) 
                        end
                        SCORE_SS_NRR(:,nrr) = SCORE_SS;
                    end
                    SCORE_SS_NRR_NRT(:,:,nrt) = SCORE_SS_NRR;
                end
                SCORE_SS_NRR_NRT_SN(:,:,:,sn) = SCORE_SS_NRR_NRT;
            end
            SCORE_SS_NRR_NRT_SN_NL(:,:,:,:,nl) = SCORE_SS_NRR_NRT_SN;
        end
        SCORE_SS_NRR_NRT_SN_NL_SQ(:,:,:,:,:,sq) = SCORE_SS_NRR_NRT_SN_NL;
    end
    SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd) = SCORE_SS_NRR_NRT_SN_NL_SQ;
end
mytime = toc;
save('result_pdf_PE.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
    'SS','NRT','NL','SN','SQ','alpha_list','sqenc_list','sensor_list');

% 
% 
%% plot Dimensionality-ROC
% SCORE_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ,2);
% figure('position',[100 100 600 300]);
% p = plot([1,5:5:300],squeeze(AUC_avg(1,1,:)));
% p.Color = 'k';
% p.LineStyle = '-';
% p.LineWidth = 2;
% p.Marker = 'o';
% p.MarkerSize = 3;
% p.MarkerEdgeColor = 'k';
% p.MarkerFaceColor = 'k';
% xlim([1,inf]);
% xticks([1,50:50:300]);
% xlabel('Dimensionality of PDF');
% ylabel('ROC-AUC');
% 
% 
%% boxplot
% SCORE_avg = squeeze(mean(SCORE_SS_NRR_NRT_SN_NL_SQ,2));
% SCORE_rep_errorbar = SCORE_avg(:,:,:);
% figure('position',[0 0 250 300])
% boxplot(squeeze(SCORE_rep_errorbar(1,1:25,:)),...
%     'symbol','','Widths',0.3);hold on;
% points = findobj(gcf, 'type', 'circle', 'Tag', 'Median');
% set(points, 'Color', 'r');
% for j=1:size(SCORE_rep_errorbar,3)
%     hold on
%     scatter(j*ones(25,1),...
%         squeeze(SCORE_rep_errorbar(1,26:50,j)),10)
% end
% xticks(1:4)
% xticklabels({'0%','15%','30%','45%'})
% % ylim([-inf,max(SCORE_rep_errorbar(:))])
% ylim([1.5,6.5])
% xlabel('Noise level')
% %     if i == 1
% ylabel('DI for Sensor 1');
% % 
% % 
% %% AUC single
% SCORE_avg = squeeze(mean(SCORE_SS_NRR_NRT_SN_NL_SQ,2));
% figure('position',[0,0,300,300]);
% labels = (1:NRT)>(NRT/2);
% i=1;
% for k = 1:size(SCORE_avg,3)% noise level
% %     for j = 1:size(SCORE_avg,3)% test
%         [FPR,TPR,~,auc] = perfcurve(labels,SCORE_avg(i,:,k),1);
% %             if auc < 0.5
% %                 auc = 1 - auc;
% %             end
%         AUC(i,k) = auc;
% %     end
%     plot(FPR,TPR,'LineWidth',2);hold on
% end
% hold off
% legend(['Noise Lev 0, AUC = ',num2str(mean(AUC(i,1),2),'%4.3f')],...
%     ['Noise Lev 1, AUC = ',num2str(mean(AUC(i,2),2),'%4.3f')],...
%     ['Noise Lev 2, AUC = ',num2str(mean(AUC(i,3),2),'%4.3f')],...
%     ['Noise Lev 3, AUC = ',num2str(mean(AUC(i,4),2),'%4.3f')],...
%     'Location','southeast','FontSize',8);
% xlabel('FPR');ylabel('TPR');
% title(['Sensor ',num2str(1)])
% % 
% AUC
SCORE_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,SN,NL,SQ,DD);
for m = 1:SS
    for i = 1:SN
        for j = 1:NL
            for k = 1:SQ
                for dd = 1:DD
                    [~,~,~,auc] = perfcurve(labels,...
                        squeeze(SCORE_avg(m,1,:,i,j,k,dd)),1);
                    AUC(m,i,j,k,dd) = auc;
                end
            end
        end
    end
end
pdfm784 = squeeze(mean(AUC(:,:,:,:,1:5),5))
pdfs784 = squeeze(std(AUC(:,:,:,:,1:5),0,5))
pdfm768 = squeeze(mean(AUC(:,:,:,:,6:10),5))
pdfs768 = squeeze(std(AUC(:,:,:,:,6:10),0,5))


% scatter(ones(25,1),SCORE_avg(1,:,1:25));hold on;
% scatter(ones(25,1),SCORE_avg(1,:,26:50));legend('1','2')