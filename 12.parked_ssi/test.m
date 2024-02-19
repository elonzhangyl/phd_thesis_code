clc
clear
%% Load data set:
dam = 0.48;
load('und_seed101To400_1e6.mat')
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

load(['dam_scour',num2str(dam),'m_seed201To300_1e6.mat'])
dataset2 = dataset2(120*20+1:end,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));

% wd = [];
% [pxx1,f1] = pwelch(dataset21(:,1,1),wd,[],[],20);
% [pxx2,f2] = pwelch(dataset21(:,1,21),wd,[],[],20);
% plot(f1,log10(pxx1),f2,log10(pxx2))
% legend('1','2')

%% main
wd_list = 10;%[0.7,0.8,0.9,1];
pca_order = 8;
sqenc_list = 130;%5:5:80;
alpha_list = 0.032;
% alpha_list = 0.0565;

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
% mu = zeros(pca_order,1);
% Sigma = eye(pca_order) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
node = 1;
resam = 1;
perid = 9;
% load('xre');
% x_re = x_re_de;
for dd = 1:DD
    dd
    % select Monto Carlo seeds
    dataset_ref = dataset1(1:resam:20*perid*60,node,...
        [1:10,201:210]+(dd-1)*10);
    dataset_test = cat(3,dataset1(1:resam:20*perid*60,...
        node,[101:110,301:310]+(dd-1)*10),...
        dataset21(1:resam:20*perid*60,node,...
        [1:10,101:110]+(dd-1)*10));
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
                parfor nrt=1:NRT % NUM RECORD TEST
                    nrt;
                    Dataset2 = dataset_test(:,nrt); %choose record 
%                     rng(100*(dd-1)+NRR+nrt);
                    rng(floor(abs(1e7*Dataset2(1))));
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset20(1,1)
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_WW_NRR = zeros(WW,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,nrr); %choose record
                         %choose record
%                         rng(100*(dd-1)+nrr);
                        rng(floor(abs(1e7*Dataset1(1))));
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
                        for ww = 1:WW % width
%                             wd = wd_list(ww);
                            SCORE_WW(ww) = KLDiv(X2,X1);%de in front
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
% save(['results_owt_KLD_',num2str(dam),...
%     '_pc',num2str(pca_order),'_sq',num2str(sqenc_list),...
%     '_de_dd10_alpha0.565.mat'],'SCORE_WW_NRR_NRT_SN_NL_SQ_DD');


% asss = squeeze(SCORE_WW_NRR_NRT_SN_NL_SQ_DD);
% NRT = 40
labels = (1:NRT)>(NRT/2);
AUC = zeros(WW,NRR,SN,NL,SQ,DD);
avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD,2));
Auc = zeros(SQ,DD);
for sq=1:SQ
for dd=1:DD
    [~,~,~,auc] = perfcurve(labels,squeeze(avg(:,dd)),1);
    Auc(sq,dd) = auc; 
end
end
mean(Auc)
