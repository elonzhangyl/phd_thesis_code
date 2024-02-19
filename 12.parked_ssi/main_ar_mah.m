clc
clear
%% Load data set:
% load('und_seed101To111.mat')
load('und_seed101To150_1e6.mat')
% dataset1=dataset1(:,1,1:20);
% load('und_seed101To111_lowpass0.1Hz_4.5e6.mat')
% load('und_seed101To111_6e5.mat')
% load('dam_scour0.12m_seed201To207.mat')
% init = 8000;
% dataset1 = dataset1;
% dataset1 = dataset1(init:init+600*20,:,:);
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

% load('dam_scour0.6m_seed201To210.mat')
load('dam_scour0.36m_seed201To225_1e6.mat')
dataset2 = dataset2(120*20+1:end,:,:);
% load('dam_scour0.36m_seed201To210_lowpass0.1Hz_4.5e6.mat')
% load('dam_scour0.6m_seed201To210_6e5.mat')
% dataset21 = dataset2;
% dataset21 = dataset2(init:init+600*20,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));


% wd = [];
% [pxx1,f1] = pwelch(dataset21(:,1,1),wd,[],[],20);
% [pxx2,f2] = pwelch(dataset21(:,1,21),wd,[],[],20);
% plot(f1,log10(pxx1),f2,log10(pxx2))
% legend('1','2')

%% sequence density ratio KL divergence for localization
node = 1;
resam = 1;
% both
% dataset_ref = dataset1(1:resam:end,node,[1:10,51:60]);
% dataset_test = cat(3,dataset1(1:resam:end,node,[11:20,61:70]),...
%     dataset21(1:resam:end,node,[1:10,21:30]));
% dataset_ref = dataset1(1:resam:end,node,[21:30,71:80]);
% dataset_test = cat(3,dataset1(1:resam:end,node,[31:40,81:90]),...
%     dataset21(1:resam:end,node,[11:20,31:40]));
% dataset_ref = dataset1(1:resam:end,node,[21:30,71:80]);
% dataset_test = cat(3,dataset1(1:resam:end,node,[31:40,81:90]),...
%     dataset21(1:resam:end,node,[1:10,11:20]));
% dataset_ref = dataset1(1:resam:end,node,1:20);
% dataset_test = cat(3,dataset1(1:resam:end,node,21:40),...
%     dataset21(1:resam:end,node,1:20));
dataset_ref = dataset1(1:resam:end,node,41:60);
dataset_test = cat(3,dataset1(1:resam:end,node,61:80),...
    dataset21(1:resam:end,node,21:40));
%2nd
% dataset_ref = dataset1(1:resam:end,node,51:70);
% dataset_test = cat(3,dataset1(1:resam:end,node,71:90),...
%     dataset21(1:resam:end,node,21:40));
% 1st
% dataset_ref = dataset1(1:resam:end,node,1:20);
% dataset_test = cat(3,dataset1(1:resam:end,node,21:40),...
%     dataset21(1:resam:end,node,1:20));
%% main
ar_order = 250;
sqenc_list = 5;%5:5:80;
alpha_list = 0;
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
                        Mah_d = mahal(squeeze(AR_Coef2(:,ss)'),...
                            squeeze(AR_Coef1(:,ss,:))');%damaged is before.
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
% figure('position',[100 100 600 150])
% plot(SCORE_SS_NRT_NS_NL_SQ_DD(:),'.','MarkerSize',20);
% % 
% % score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
% % plot(score_avg(:),'.');
% 
% 
% labels = (1:NRT)>10;
% 
% [~,~,~,auc] = perfcurve(labels,...
%     (SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:)),1)


% NRT = 20
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NS,NL,SQ,DD);

for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for ss = 1:SS
                    [~,~,~,auc] = perfcurve(labels,...
                        squeeze(SCORE_SS_NRT_NS_NL_SQ_DD...
                        (ss,:,ns,nl,sq,dd)),1);
                    AUC(ss,ns,nl,sq,dd) = auc;
                end
            end
        end
    end
end 
squeeze(AUC)



