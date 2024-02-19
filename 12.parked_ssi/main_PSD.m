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
load('dam_scour0.24m_seed201To225_1e6.mat')
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
dataset_ref = dataset1(1:resam:end,node,1:20);
dataset_test = cat(3,dataset1(1:resam:end,node,21:40),...
    dataset21(1:resam:end,node,1:20));
% dataset_ref = dataset1(1:resam:end,node,41:60);
% dataset_test = cat(3,dataset1(1:resam:end,node,61:80),...
%     dataset21(1:resam:end,node,21:40));
%2nd
% dataset_ref = dataset1(1:resam:end,node,51:70);
% dataset_test = cat(3,dataset1(1:resam:end,node,71:90),...
%     dataset21(1:resam:end,node,21:40));
% 1st
% dataset_ref = dataset1(1:resam:end,node,1:20);
% dataset_test = cat(3,dataset1(1:resam:end,node,21:40),...
%     dataset21(1:resam:end,node,1:20));

%% main
wd_list = 0.7;%[0.7,0.8,0.9,1];
sqenc_list = 20;%5:20;
alpha_list = 0.2;
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;%seed
NRR = size(dataset_ref,2);
NRT = size(dataset_test,2); 
DD = length(wd_list);
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
tic
% mu = zeros(sqenc_list,1);
% Sigma = eye(sqenc_list) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
for dd = 1:DD
    wd = wd_list(dd);
    SCORE_SS_NRR_NRT_SN_NL_SQ = zeros(SS,NRR,NRT,SN,NL,SQ);
    for sq = 1:SQ
        sqenc = sqenc_list(sq);
        SCORE_SS_NRR_NRT_SN_NL = zeros(SS,NRR,NRT,SN,NL);
        for nl = 1:NL % noise level
%             nla
            alpha =alpha_list(nl);
            SCORE_SS_NRR_NRT_SN = zeros(SS,NRR,NRT,SN);
            for sn = 1:SN % seed noise
                SCORE_SS_NRR_NRT = zeros(SS,NRR,NRT);
                for nrt=1:NRT % NUM RECORD TEST
                    nrt
                    Dataset2 = dataset_test(:,nrt); %choose record 
    %                 rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_SS_NRR = zeros(SS,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,nrr); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = zscore(Dataset10); %normlization
                        SCORE_SS = zeros(SS,1);
                        for ss = 1:SS
                            [pxx1,fxx1] = pwelch(Dataset100,[],[],[]);
                            [pxx2,fxx2] = pwelch(Dataset200,[],[],[]);
%                             plot(fxx1,log10(pxx1),fxx2,log10(pxx2))
                            SCORE_SS(ss) = mean(pxx1./pxx2)
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
% save('results_10m_sq80.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
%     'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');
% 
% figure('position',[100 100 600 150])
% plot(SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:),'.','MarkerSize',20);
% % 
% % score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
% % plot(score_avg(:),'.');
% 
% 
% labels = (1:NRT)>10;
% 
% [~,~,~,auc] = perfcurve(labels,...
%     (SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:)),1)


% labels = (1:NRT)>(NRT/2);
% AUC = zeros(SS,NRR,SN,NL,SQ,DD);
% 
% for dd = 1:DD
%     for sq = 1:SQ
%         for nl = 1:NL % noise level
%             for sn = 1:SN % seed noise
%                 for nrr = 1:NRR %NUM RECORD reference
%                     for ss = 1:SS
%                         [~,~,~,auc] = perfcurve(labels,...
%                             squeeze(SCORE_SS_NRR_NRT_SN_NL_SQ_DD...
%                             (ss,nrr,[1:10,11:20],sn,nl,sq,dd)),1);
%                         AUC(ss,nrr,sn,nl,sq,dd) = auc;
%                     end
%                 end
%             end
%         end
%     end
% end 
% squeeze(AUC)

labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NRR,SN,NL,SQ,DD);
score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);

for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for nrr = 1 %NUM RECORD reference
                    for ss = 1:SS
                        [~,~,~,auc] = perfcurve(labels,...
                            squeeze(score_avg),1);
                        AUC(ss,nrr,sn,nl,sq,dd) = auc;
                    end
                end
            end
        end
    end
end 
squeeze(AUC)




