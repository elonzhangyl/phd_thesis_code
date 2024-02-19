clc
clear
%% Load data set:
% load('und_seed101To111.mat')
load('und_seed101To150_1e6.mat')
% load('und_seed101To200_1e6.mat')
% dataset1=dataset1(:,1,1:20);
% load('und_seed101To111_lowpaww0.1Hz_4.5e6.mat')
% load('und_seed101To111_6e5.mat')
% load('dam_scour0.12m_seed201To207.mat')
% init = 8000;
% dataset1 = dataset1;
% dataset1 = dataset1(init:init+600*20,:,:);
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

% load('dam_scour0.6m_seed201To210.mat')
load('dam_scour0.36m_seed201To225_1e6.mat')
% load('dam_scour0.36m_seed201To250_1e6.mat')
dataset2 = dataset2(120*20+1:end,:,:);
% load('dam_scour0.36m_seed201To210_lowpaww0.1Hz_4.5e6.mat')
% load('dam_scour0.6m_seed201To210_6e5.mat')
% dataset21 = dataset2;
% dataset21 = dataset2(init:init+600*20,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));

% wd = [];
% [pxx1,f1] = pwelch(dataset21(:,1,1),wd,[],[],20);
% [pxx2,f2] = pwelch(dataset21(:,1,21),wd,[],[],20);
% plot(f1,log10(pxx1),f2,log10(pxx2))
% legend('1','2')

%% main
wd_list = 10;%[0.7,0.8,0.9,1];
pca_order = 8;
sqenc_list = 180;%5:5:80;
alpha_list = 0;
SQ = length(sqenc_list);
NL = length(alpha_list);
SN = 1;%seed
NRR = 20; % number of reference time series
NRT = 40; % number of test time series
WW = length(wd_list);
DD = 1;
SCORE_WW_NRR_NRT_SN_NL_SQ_DD = zeros(WW,NRR,NRT,SN,NL,SQ,DD);
tic
% mu = zeros(sqenc_list,1);
% Sigma = eye(sqenc_list) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
for dd = 1:DD
    dd
    % select Monto Carlo seeds
    node = 1;
    resam = 1;
%     dataset_ref = dataset1(1:resam:end,node,[1:10,101:110]+(dd-1)*10);
%     dataset_test = cat(3,dataset1(1:resam:end,node,[51:60,151:160]+(dd-1)*10),...
%         dataset21(1:resam:end,node,[1:10,51:60]+(dd-1)*10));
dataset_ref = dataset1(1:resam:end,node,[1:10,51:60]);
dataset_test = cat(3,dataset1(1:resam:end,node,[11:20,61:70]),...
    dataset21(1:resam:end,node,[1:10,21:30]));
    % run
    SCORE_WW_NRR_NRT_SN_NL_SQ = zeros(WW,NRR,NRT,SN,NL,SQ);
    for sq = 1:SQ
        sqenc = sqenc_list(sq);
        SCORE_WW_NRR_NRT_SN_NL = zeros(WW,NRR,NRT,SN,NL);
        for nl = 1:NL % noise level
%             nla
            alpha =alpha_list(nl);
            SCORE_WW_NRR_NRT_SN = zeros(WW,NRR,NRT,SN);
            for sn = 1:SN % seed noise
                SCORE_WW_NRR_NRT = zeros(WW,NRR,NRT);
                for nrt=1:NRT % NUM RECORD TEST
                    nrt
                    Dataset2 = dataset_test(:,:,nrt); %choose record 
    %                 rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_WW_NRR = zeros(WW,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,:,nrr); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = zscore(Dataset10); %normlization
                        SCORE_WW = zeros(WW,1);

                        X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                        X2 = zeros(size(Dataset200,1)-sqenc,sqenc);
                        for j=1:(size(Dataset100,1)-sqenc)
                            X1(j,:) = Dataset100(j:j+sqenc-1);
                        end
                        for j=1:(size(Dataset200,1)-sqenc)
                            X2(j,:) = Dataset200(j:j+sqenc-1);
                        end
                        [coeff,score,latent] = pca(X1);
                        X1_PCA = score(:,1:pca_order);
                        X2_PCA = X2 * coeff(:,1:pca_order);
%                             else
%                             X2_PCA = X2 * coeff(:,1:pca_order);
%                             end
                        X1_PCA = zscore(X1_PCA);
                        X2_PCA = zscore(X2_PCA);
%                             [~,score1,~] = pca(X1);
%                             [~,score2,~] = pca(X2);
%                             X1_PCA = score1(:,1:pca_order);
%                             X2_PCA = score2(:,1:pca_order);  
                        for ww = 1:WW % width
                            wd = wd_list(ww);
                            x_de=X1_PCA';
                            x_nu=X2_PCA';
%                             uLSIF
                            % wd=10;lambda=0.01;b=500;
                            [~,wh_x_re1,~,~] = KLIEP2(x_de,x_nu,x_nu,wd); 
%                             [~,wh_x_re1] = uLSIF(x_de,x_nu,x_nu,wd,lambda,b,5);
                            SCORE_WW(ww) = mean(log(wh_x_re1))
%                             SCORE_WW(ww) = KLDiv(X1_PCA,X2_PCA)+KLDiv(X2_PCA,X1_PCA)
                            x_de=X2_PCA';
                            x_nu=X1_PCA';
                            [~,wh_x_re2,~,~] = KLIEP2(x_de,x_nu,x_nu,wd);
%                             [~,wh_x_re2] = uLSIF(x_de,x_nu,x_nu,wd,lambda,b,5);
                            SCORE_WW(ww) = mean(log(wh_x_re2))
%                             SCORE_WW(ww) = mean(log(wh_x_re1))+mean(log(wh_x_re2))
%                             SCORE_WW(ww) = -1/2*mean((wh_x_re1).^2)+mean((wh_x_re1))+...
%                                 -1/2*mean((wh_x_re2).^2)+mean((wh_x_re2))
%                             SCORE_WW(ww) = -1/2*mean((wh_x_re2).^2)+mean((wh_x_re2))
%                             SCORE_WW(ww) = KLDiv(X1,X2)
                        end
                        SCORE_WW_NRR(:,nrr) = SCORE_WW;
                    end
                    SCORE_WW_NRR_NRT(:,:,nrt) = SCORE_WW_NRR;
                end
                SCORE_WW_NRR_NRT_SN(:,:,:,sn) = SCORE_WW_NRR_NRT;
            end
            SCORE_WW_NRR_NRT_SN_NL(:,:,:,:,nl) = SCORE_WW_NRR_NRT_SN;
        end
        SCORE_WW_NRR_NRT_SN_NL_SQ(:,:,:,:,:,sq) = SCORE_WW_NRR_NRT_SN_NL;
    end
    SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,dd) = SCORE_WW_NRR_NRT_SN_NL_SQ;
end
mytime = toc;
% save('results_noZscore2PCA_PCA_DRE.mat','SCORE_WW_NRR_NRT_SN_NL_SQ_DD');

% figure('position',[100 100 600 150])
% plot(SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:),'.','MarkerSize',20);
% % 
% % score_avg = mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD,2);
% % plot(score_avg(:),'.');
% 
% 
% labels = (1:NRT)>10;
% 
% [~,~,~,auc] = perfcurve(labels,...
%     (SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:)),1)


% NRT = 40
labels = (1:NRT)>(NRT/2);
AUC = zeros(WW,NRR,SN,NL,SQ,DD);
avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD,2));
Auc = zeros(DD,1);
for dd=1:DD
    [~,~,~,auc] = perfcurve(labels,squeeze(avg(:,dd)),1)
    Auc(dd) = auc; 
end
mean(Auc)
std(Auc)

% plot(squeeze(avg),'.-')
for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for nrr = 1:NRR %NUM RECORD reference
                    for ww = 1:WW
                        [~,~,~,auc] = perfcurve(labels,...
                            squeeze(SCORE_WW_NRR_NRT_SN_NL_SQ_DD...
                            (ww,nrr,:,sn,nl,sq,dd)),1);
                        AUC(ww,nrr,sn,nl,sq,dd) = auc;
                    end
                end
            end
        end
    end
end 
squeeze(AUC)



