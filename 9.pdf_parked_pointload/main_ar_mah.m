clear;clc;

tic;

%% Tuning parameter
% ar_order = 25;
% sqenc_list = 13%2:25;%5
% SQ = length(sqenc_list);
% alpha_list = [0 0.15 0.30 0.45];
% sensor_list=1:4;
% NL = length(alpha_list);
% SS = length(sensor_list);
% NS = 1;
% NRR = 25;
% NRT = 50;  
% DD = 10;

ar_order = 25;
sqenc_list = 13%2:25;%5
SQ = length(sqenc_list);
alpha_list = [0];
sensor_list= 1:4;
NL = length(alpha_list);
SS = length(sensor_list);
NS = 1;
NRR = 25;
NRT = 1;  
DD = 1;

SCORE_SS_NRT_NS_NL_SQ_DD = zeros(SS,NRT,NS,NL,SQ,DD);
for dd = 1:DD
    dd
    file_ref = {'und_seed201To205.mat',...
        'und_seed211To215.mat',...
        'und_seed221To225.mat',...
        'und_seed231To235.mat',...
        'und_seed241To245.mat',...
        'und_seed201To205.mat',...
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
        'und_seed246To250.mat',...
        'und_seed206To210.mat',...
        'und_seed216To220.mat',...
        'und_seed226To230.mat',...
        'und_seed236To240.mat',...
        'und_seed246To250.mat'};
    load(string(file_test_1(dd)));
    dataset_test_1 = dataset1;
    dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);


    file_test_2 = {'dam_0.02_seed301To305.mat',...
        'dam_0.02_seed306To310.mat',...
        'dam_0.02_seed311To315.mat',...
        'dam_0.02_seed316To320.mat',...
        'dam_0.02_seed321To325.mat',...
        'dam_0.04_seed401To405.mat',...
        'dam_0.04_seed406To410.mat',...
        'dam_0.04_seed411To415.mat',...
        'dam_0.04_seed416To420.mat',...
        'dam_0.04_seed421To425.mat'};
    load(string(file_test_2(dd)));

    dataset_test_2 = dataset2;
    dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);
    dataset_test = cat(2,dataset_test_1,dataset_test_2);

  
    %% Damage Index
    SCORE_SS_NRT_NS_NL_SQ = zeros(SS,NRT,NS,NL,SQ);
    for sq = 1:SQ
        sq
        sqenc = sqenc_list(sq);
        CoefSel = 2:sqenc;
        SCORE_SS_NRT_NS_NL = zeros(SS,NRT,NS,NL);
        for nl = 1:NL % noise level
%         nl
        alpha =alpha_list(nl);
        SCORE_SS_NRT_NS = zeros(SS,NRT,NS);
            for ns = 1:NS %  noise seed tier                        
                AR_Coef1 = zeros(length(CoefSel),NRR,SS);
                for nrr = 1:NRR %num record referece
                    Dataset1 = dataset_ref(:,nrr,:);% choose record
                    rng(floor(abs(1e7*Dataset1(1))));
                    Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                        .*(ones(size(Dataset1)).*rms(Dataset1));%add noise
                    Dataset100 = zscore(Dataset10);%normlization
                    %% AR Coefficients
                    for ss = 1:SS%sensors
                        m1 = ar(Dataset100(:,ss),ar_order,'burg');
                        AR_Coef1(:,nrr,ss) = m1.A(CoefSel);
                    end
                end
                %dam AR coefs
                SCORE_SS_NRT = zeros(SS,NRT);
                for nrt = 1:NRT %num of records test
%                     nrt
                    Dataset2 = dataset_test(:,nrt,:);% choose record
                    rng(floor(abs(1e7*Dataset2(1))));
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2));%add noise
                    Dataset200 = zscore(Dataset20);%normlization
                    AR_Coef2 = zeros(length(CoefSel),SS);
                    %% AR Coefficients
                    for ss = 1:SS % sensors
                        m2 = ar(Dataset200(:,ss),ar_order,'burg');
                        AR_Coef2(:,ss) = m2.A(CoefSel);
                    end
                    %% Mahalanobis Distance
                    SCORE_SS = zeros(SS,1);
                    for ss = 1:SS % sensor number
                        Mah_d = mahal(squeeze(AR_Coef2(:,ss)'),...
                            squeeze(AR_Coef1(:,:,ss))');%damaged is before.
                        SCORE_SS(ss) = Mah_d;
                    end
                    SCORE_SS_NRT(:,nrt) = SCORE_SS;
                end
                SCORE_SS_NRT_NS(:,:,ns) = SCORE_SS_NRT;
            end
            SCORE_SS_NRT_NS_NL(:,:,:,nl) = SCORE_SS_NRT_NS;
        end
        SCORE_SS_NRT_NS_NL_SQ(:,:,:,:,sq) = SCORE_SS_NRT_NS_NL;
    end
    SCORE_SS_NRT_NS_NL_SQ_DD(:,:,:,:,:,dd) = SCORE_SS_NRT_NS_NL_SQ;
end
% mytime = toc;
% save('result_ar_mah_arOrder25_ARslct13.mat','SCORE_SS_NRT_NS_NL_SQ_DD',...
%     'SS','NRT','NL','NS','alpha_list','sensor_list','sqenc_list');
% 
% 
% 
% %% AUC
DD=10;SS=4;NS=1;NL=4;SQ=1;DD=10;
SCORE_avg = (SCORE_SS_NRT_NS_NL_SQ_DD(:,:,:,:,:,:));
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NS,NL,SQ,DD);
for m = 1:SS
    for i = 1:NS
        for j = 1:NL
            for sq = 1:SQ
                for dd = 1:DD
                    [~,~,~,auc] = perfcurve(labels,...
                        squeeze(SCORE_avg(m,:,i,j,sq,dd)),1);
                    AUC(m,i,j,sq,dd) = auc;
                end
            end
        end
    end
end
% arm784 = squeeze(mean(AUC(:,:,:,1,1:5),5))
% ars784 = squeeze(std(AUC(:,:,:,1,1:5),0,5))
% arm768 = squeeze(mean(AUC(:,:,:,1,6:10),5))
% ars768 = squeeze(std(AUC(:,:,:,1,6:10),0,5))
% 
% % best k=6
% % aa = mean(AUC(1,1,1,:,:),5);
% % plot(1:length(aa),squeeze(aa))
% 
% 



