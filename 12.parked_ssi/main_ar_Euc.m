clc
clear
%% Load data set:
load('und_seed101To111.mat')
% load('und_seed101To111_6e5.mat')
% load('dam_scour0.12m_seed201To207.mat')
dataset1 = dataset1(6000:6000+100*20,:,:);
load('dam_scour0.6m_seed201To210.mat')
% dataset21 = dataset2;
dataset21 = dataset2(6000:6000+100*20,:,:);
% load('dam_scour1.2m_seed201To210.mat')
% dataset22 = dataset2;
% load('dam_scour1.8m_seed201To210.mat')
% dataset23 = dataset2;
% load('dam_scour2.4m_seed201To210.mat')
% dataset24 = dataset2;
% load('dam_scour3.6m_seed201To210.mat')
% dataset25 = dataset2;
% load('dam_scour4.8m_seed201To210.mat')
% dataset26 = dataset2;

%% sequence density ratio KL divergence for localization
node = 1;
% dataset11 = zeros(size(dataset1));
% for i = 1:10
% dataset11(:,node,i) = dataset1(:,node,i) + 0.3*randn(size(dataset1(:,node,i)))...
%                         .*(ones(size(dataset1(:,node,i))).*rms(dataset1(:,node,i),1));
% std(dataset11(:,node,i));
% end
% dataset1 = dataset11;

resam = 1;
dataset_ref = dataset1(1:resam:end,node,1);
% dataset_test = cat(3,dataset1(:,node,11:20),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10),...
%     dataset25(:,node,1:10),dataset26(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,2:11),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10));
dataset_test = cat(3,dataset1(1:resam:end,node,2:11),dataset21(1:resam:end,node,1:10));

% dataset_ref = dataset_ref(1:end/2,:,:);
% dataset_test = dataset_test(1:end/2,:,:);
%% main
ar_order = 20;
sqenc_list = 7;%5:5:80;
alpha_list = 0.2;
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;%seed
NRR = size(dataset_ref,3);
NRT = size(dataset_test,3); 
DD = 1;
NS = 1;
SCORE_SS_NRT_NS_NL_SQ_DD = zeros(SS,NRT,NS,NL,SQ,DD);
tic
% mu = zeros(sqenc_list,1);
% Sigma = eye(sqenc_list) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
for dd = 1:DD
    SCORE_SS_NRT_NS_NL_SQ = zeros(SS,NRT,NS,NL,SQ);
    for sq = 1:SQ
%         sq
        sqenc = sqenc_list(sq);
        CoefSel = 2:sqenc;
        SCORE_SS_NRT_NS_NL = zeros(SS,NRT,NS,NL);
        for nl = 1:NL % noise level
%         nl
        alpha =alpha_list(nl);
        SCORE_SS_NRT_NS = zeros(SS,NRT,NS);
            for ns = 1:NS %  noise seed tier                        
                AR_Coef1 = zeros(length(CoefSel),SS,NRR);
                for nrr = 1:NRR %num record referece
                    Dataset1 = dataset_ref(:,:,nrr);% choose record
    %                 rng(100*(ns-1)+nrr)
                    Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                        .*(ones(size(Dataset1)).*rms(Dataset1));%add noise
                    Dataset100 = zscore(Dataset10);%normlization
                    %% AR Coefficients
                    for ss = 1:SS%sensors
                        m1 = ar(Dataset100(:,ss),ar_order,'burg');
                        AR_Coef1(:,ss,nrr) = m1.A(CoefSel);
                    end
                end
                %dam AR coefs
                SCORE_SS_NRT = zeros(SS,NRT);
                for nrt = 1:NRT %num of records test
                    nrt
                    Dataset2 = dataset_test(:,:,nrt);% choose record
    %                 rng(100*(ns-1)+NRR+nrt)
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
                        Mah_d = norm(squeeze(AR_Coef2(:,ss)')-...
                            squeeze(AR_Coef1(:,ss,:))');%damaged is before                        
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
mytime = toc;
% save('results_10m_sq80.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
%     'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');
% 
figure('position',[100 100 600 150])
plot(SCORE_SS_NRT_NS_NL_SQ_DD(:),'.','MarkerSize',20);
% % 
% % score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
% % plot(score_avg(:),'.');
% 
% 
% labels = (1:NRT)>10;
% 
% [~,~,~,auc] = perfcurve(labels,...
%     (SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:)),1)


NRT = 20
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NS,NL,SQ,DD);

for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for ss = 1:SS
                    [~,~,~,auc] = perfcurve(labels,...
                        squeeze(SCORE_SS_NRT_NS_NL_SQ_DD...
                        (ss,[1:10,11:20],ns,nl,sq,dd)),1);
                    AUC(ss,ns,nl,sq,dd) = auc;
                end
            end
        end
    end
end 
squeeze(AUC)



