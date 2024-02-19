clc
clear
%% Load data set:
load('und_seed101To125.mat')
% load('dam_scour0.12m_seed201To207.mat')
load('dam_scour0.6m_seed201To210.mat')
dataset21 = dataset2;
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

% load('dam_scour5m_seed201To225.mat')
% load('dam_scour7.8m_seed201To210.mat')
% load('dam_scour10m_seed201To210.mat')

% 
% wd = 2048;
% [pxx1,f1] = pwelch(dataset1(:,5,1),wd,[],[],20);
% [pxx2,f2] = pwelch(dataset2(:,5,1),wd,[],[],20);
% 
% 
% plot(f1,log10(pxx1),f2,log10(pxx2))
% legend('1','2')
% 
% spectrogram(dataset(:,5,93),[],[],[],320,'yaxis')
% 
% qqplot(dataset1(1:end/10,6,2))
% for i = 1:6
%     x = dataset2(:,6,i);
%     xx = x + 0*randn(size(x))...
%                         .*(ones(size(x)).*rms(x,1));
%     [pxx,f] = pwelch(xx,[],[],[],20);
%     freq(i) = f(find(pxx == max(pxx(1:floor(end/10)))));
% end
% mean(freq)
%% sequence density ratio KL divergence for localization
node = 1;
% dataset11 = zeros(size(dataset1));
% for i = 1:10
% dataset11(:,node,i) = dataset1(:,node,i) + 0.3*randn(size(dataset1(:,node,i)))...
%                         .*(ones(size(dataset1(:,node,i))).*rms(dataset1(:,node,i),1));
% std(dataset11(:,node,i));
% end
% dataset1 = dataset11;

dataset_ref = dataset1(:,node,11);
% dataset_test = cat(3,dataset1(:,node,11:20),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10),...
%     dataset25(:,node,1:10),dataset26(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,2:11),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10));
dataset_test = cat(3,dataset1(:,node,16:25),dataset21(:,node,1:10));
dataset_ref = dataset_ref(1:end/2,:,:);
dataset_test = dataset_test(1:end/2,:,:);

%% main
wd_list = 0.8;%[0.7,0.8,0.9,1];
pca_order = 20;
sqenc_list = 160;%5:5:80;
alpha_list = 0;
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;%seed
NRR = size(dataset_ref,3);
NRT = size(dataset_test,3); 
DD = length(wd_list);
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
tic
% mu = zeros(sqenc_list,1);
% Sigma = eye(sqenc_list) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
%% pca
for dd = 1:DD
    for sq = 1:SQ
        sqenc = sqenc_list(sq);
        for nl = 1:NL % noise level
            alpha =alpha_list(nl);
            for sn = 1:SN % seed noise
                for nrt=1:NRT % NUM RECORD TEST
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,:,nrr); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = zscore(Dataset10); %normlization
                        SCORE_SS = zeros(SS,1);
                        for ss = 1:SS
                            X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                            for j=1:(size(Dataset100,1)-sqenc)
                                X1(j,:) = Dataset100(j:j+sqenc-1,ss);
                            end
                                [coeff,score,latent] = pca(X1);
                                X1_PCA = score(:,1:pca_order);
                                X1_PCA = zscore(X1_PCA);
                        end
                    end
                end
            end
        end
    end
end
%% main
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
                parfor nrt=1:NRT % NUM RECORD TEST
                    nrt
                    Dataset2 = dataset_test(:,:,nrt); %choose record 
    %                 rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_SS_NRR = zeros(SS,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,:,nrr); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = zscore(Dataset10); %normlization
                        SCORE_SS = zeros(SS,1);
                        for ss = 1:SS
                            X2 = zeros(size(Dataset200,1)-sqenc,sqenc);
                            for j=1:(size(Dataset200,1)-sqenc)
                                X2(j,:) = Dataset200(j:j+sqenc-1,ss);
                            end
                            
                            X2_PCA = X2 * coeff(:,1:pca_order);
                            X2_PCA = zscore(X2_PCA);
%                             [~,score1,~] = pca(X1);
%                             [~,score2,~] = pca(X2);
%                             X1_PCA = score1(:,1:pca_order);
%                             X2_PCA = score2(:,1:pca_order);  
                            x_de=X1_PCA';
                            x_nu=X2_PCA';
                            [~,wh_x_re,~,~] = KLIEP2(x_de,x_nu,x_nu,wd);              
                            SCORE_SS(ss) = mean(log(wh_x_re))
                            
%                             SCORE_SS(ss) = KLDiv(X1,X2)
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
figure('position',[100 100 600 150])
plot(SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:),'.','MarkerSize',20);
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
AUC = zeros(SS,NRR,SN,NL,SQ,DD);

for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for nrr = 1:NRR %NUM RECORD reference
                    for ss = 1:SS
                        [~,~,~,auc] = perfcurve(labels,...
                            squeeze(SCORE_SS_NRR_NRT_SN_NL_SQ_DD...
                            (ss,nrr,[1:10,11:20],sn,nl,sq,dd)),1);
                        AUC(ss,nrr,sn,nl,sq,dd) = auc;
                    end
                end
            end
        end
    end
end 
squeeze(AUC)



