clear;clc;

tic;
%% Tuning parameter
sqenc_list = [1 5:5:250];%[1 5:5:300];
alpha_list = [0];
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;
NRR = 25;
NRT = 50; 
DD = 5;
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
for dd = 1:DD
    dd
    file_ref = {'und_seed201To205.mat',...
        'und_seed211To215.mat',...
        'und_seed221To225.mat',...
        'und_seed231To235.mat',...
        'und_seed241To245.mat'};
    load(string(file_ref(dd)));
    dataset_ref=dataset1;
    dataset_ref = reshape(dataset_ref,1200*20,25,4);
    
    file_test_1 = {'und_seed206To210.mat',...
        'und_seed216To220.mat',...
        'und_seed226To230.mat',...
        'und_seed236To240.mat',...
        'und_seed246To250.mat'};
    load(string(file_test_1(dd)));
    dataset_test_1 = dataset1;
    dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);


    file_test_2 = {'dam_0.01_seed1001To1005.mat',...
        'dam_0.01_seed1006To1010.mat',...
        'dam_0.01_seed1011To1015.mat',...
        'dam_0.01_seed1016To1020.mat',...
        'dam_0.01_seed1021To1025.mat'
%         'dam_0.02_seed301To305.mat',...
%         'dam_0.02_seed306To310.mat',...
%         'dam_0.02_seed311To315.mat',...
%         'dam_0.02_seed316To320.mat',...
%         'dam_0.02_seed321To325.mat'
%         'dam_0.03_seed3001To3005.mat',...
%         'dam_0.03_seed3006To3010.mat',...
%         'dam_0.03_seed3011To3015.mat',...
%         'dam_0.03_seed3016To3020.mat',...
%         'dam_0.03_seed3021To3025.mat'
%         'dam_0.04_seed401To405.mat',...
%         'dam_0.04_seed406To410.mat',...
%         'dam_0.04_seed411To415.mat',...
%         'dam_0.04_seed416To420.mat',...
%         'dam_0.04_seed421To425.mat'
        };
    load(string(file_test_2(dd)));

    dataset_test_2 = dataset2;
    dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);
    dataset_test = cat(2,dataset_test_1,dataset_test_2);

  
    %% Damage Index
    SCORE_SS_NRR_NRT_SN_NL_SQ = zeros(SS,NRR,NRT,SN,NL,SQ);
    parfor sq = 1:SQ
        sq
        sqenc = sqenc_list(sq);
        SCORE_SS_NRR_NRT_SN_NL = zeros(SS,NRR,NRT,SN,NL);
        for nl = 1:NL % noise level
%             nl
            alpha =alpha_list(nl);
            SCORE_SS_NRR_NRT_SN = zeros(SS,NRR,NRT,SN);
            for sn = 1:SN % seed noise
%                 sn
                SCORE_SS_NRR_NRT = zeros(SS,NRR,NRT);
                for nrt=1:NRT % NUM RECORD TEST
%                     nrt
                    Dataset2 = dataset_test(:,nrt,:); %choose record 
                    rng(floor(abs(1e7*Dataset2(1))));
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_SS_NRR = zeros(SS,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,nrr,:); %choose record
                         %choose record
                        rng(floor(abs(1e7*Dataset1(1))));
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
                        SCORE_SS(ss) = KLDiv(X1,X2);
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
% save('result_pdf_KLDiv_tuning_0.01.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
%     'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');

DD=1;
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


%% plot Dimensionality-ROC
figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.5]*29.1);
% axes('position',[0.1 0.1 0.9 0.9]);
xx = [1 5:5:250];
yy1 = squeeze(mean(AUC(:,:,:,:,16:20),5));
yy2 = squeeze(mean(AUC(:,:,:,:,11:15),5));
yy3 = squeeze(mean(AUC(:,:,:,:,6:10),5));
yy4 = squeeze(mean(AUC(:,:,:,:,1:5),5));
p1 = plot(xx,yy1);hold on
% p1.Color = [0, 0.4470, 0.7410];
p1.LineStyle = '-';
p1.LineWidth = 1;
p1.Marker = 'o';
p1.MarkerSize = 2.5;
% p1.MarkerEdgeColor = [0, 0.4470, 0.7410];
p1.MarkerFaceColor = 'w';
p2 = plot(xx,yy2);hold on
% p2.Color = [0.8500, 0.3250, 0.0980];
p2.LineStyle = '-';
p2.LineWidth = 1;
p2.Marker = 'v';
p2.MarkerSize = 2.5;
% p2.MarkerEdgeColor = [0.8500, 0.3250, 0.0980];
p2.MarkerFaceColor = 'w';
p3 = plot(xx,yy3);hold on
% p3.Color = [0, 0.4470, 0.7410];
p3.LineStyle = '-';
p3.LineWidth = 1;
p3.Marker = '*';
p3.MarkerSize = 2.5;
% p3.MarkerEdgeColor = [0, 0.4470, 0.7410];
p3.MarkerFaceColor = 'w';
p4 = plot(xx,yy4);hold off
% p4.Color = [0.8500, 0.3250, 0.0980];
p4.LineStyle = '-';
p4.LineWidth = 1;
p4.Marker = '+';
p4.MarkerSize = 2.5;
% p4.MarkerEdgeColor = [0.8500, 0.3250, 0.0980];
p4.MarkerFaceColor = 'w';
xlim([1,200]);
% ylim([0 1.05]);
xticks([1,50:50:300]);
xlabel('Dimensionality of PDF, {\it{\fontname{Times} k}}');
ylabel('Average ROC-AUC');
% ytickformat('percentage')
% yline(0.05,'r');
% xlabel('Dimensionality of PDF, $k$','Interpreter', 'latex');
% ylabel('ROC-AUC');
legend('4% stiffness reduction',...
    '3% stiffness reduction',...
    '2% stiffness reduction',...
    '1% stiffness reduction',...
    'location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig4.Dimensionality-ROC.eps','Resolution',1000)

%% plot
score_tuning = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
SCORE_avg = squeeze(score_tuning);
SKL = zeros(SQ,1);
for sq = 1:SQ
    skl = SCORE_avg(:,sq,:);
%         std(skl)
    th = prctile(skl(:),95)
    SKL(sq) = sum(skl(:) > th); 
end
plot(1:length(SKL),SKL)
ylim([0 0.1]);
for i = 1:SQ
    subplot(6,6,i)
    aa = score_tuning(:,:,:,:,:,i,:);
    histogram(aa(:),20)
end

Y = prctile(X,p1)