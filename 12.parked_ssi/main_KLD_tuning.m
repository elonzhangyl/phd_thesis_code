clc
clear
%% Load data set:
% load('und_seed101To150_1e6.mat')
load('und_seed101To320_1e6.mat')
% dataset1=dataset1(:,1,1:20);
% load('und_seed101To111_lowpaww0.1Hz_4.5e6.mat')
% load('und_seed101To111_6e5.mat')
% load('dam_scour0.12m_seed201To207.mat')
% init = 8000;
% dataset1 = dataset1;
% dataset1 = dataset1(init:init+600*20,:,:);
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

% load('dam_scour0.36m_seed201To250_1e6.mat')
load('dam_scour0.36m_seed201To260_1e6.mat')
% load('dam_scour1.2m_seed201To210.mat')
dataset2 = dataset2(120*20+1:end,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));

% wd = [];
% [pxx1,f1] = pwelch(dataset21(:,1,1),wd,[],[],20);
% [pxx2,f2] = pwelch(dataset21(:,1,21),wd,[],[],20);
% plot(f1,log10(pxx1),f2,log10(pxx2))
% legend('1','2')

%% main
tic
wd_list = 10;%[0.7,0.8,0.9,1];
pca_order = 8;
sqenc_list = [1:9,10:10:200];%5:5:80;
alpha_list = 0.08;
SQ = length(sqenc_list);
NL = length(alpha_list);
SN = 1;%seed
NRR = 20; % number of reference time series
NRT = 40; % number of test time series
WW = length(wd_list);
DD = 6;
SCORE_WW_NRR_NRT_SN_NL_SQ_DD = zeros(WW,NRR,NRT,SN,NL,SQ,DD);
node = 1;
resam = 1;
perid = 9;
for dd = 2
    dd
    % select Monto Carlo seeds
    if dd == 6
    dataset_ref = dataset1(1:resam:20*perid*60,node,...
        [1:10,121:130]+(dd-1)*10);
    dataset_test = cat(3,dataset1(1:resam:20*perid*60,...
        node,[61:70,181:190]+(dd-1)*10),...
        dataset21(1:resam:20*perid*60,node,...
        [1:10,61:70]+(dd-1)*10));
    else
    dataset_ref = dataset1(1:resam:20*perid*60,node,...
        [1:10,121:130]+(dd-1)*10);
    dataset_test = cat(3,dataset1(1:resam:20*perid*60,...
        node,[51:60,171:180]+(dd-1)*10),...
        dataset21(1:resam:20*perid*60,node,...
        [1:10,61:70]+(dd-1)*10));
    end
    % run
    SCORE_WW_NRR_NRT_SN_NL_SQ = zeros(WW,NRR,NRT,SN,NL,SQ);
    for sq = 1:SQ
        sq
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
                    Dataset2 = dataset_test(:,:,nrt); %choose record 
                    rng(100*(dd-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = zscore(Dataset20); %normlization 
                    SCORE_WW_NRR = zeros(WW,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,:,nrr); %choose record
                         %choose record
                        rng(100*(dd-1)+nrr);
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
                            wd = wd_list(ww);
                            SCORE_WW(ww) = KLDiv(X2,X1)%de in front
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
save('results_KLD_tuning_sq_0.36_alpha0.0565_X2X1.mat','SCORE_WW_NRR_NRT_SN_NL_SQ_DD','sqenc_list');
figure()
sqenc_list = [1:9,10:10:300];%5:5:80;
WW=1;NRR=20;NRT=40;SN=1;NL=1;SQ=length(sqenc_list);DD=6;
avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,2),2));
labels = (1:NRT)>(NRT/2);
for sq=1:SQ
    [~,~,~,auc] = perfcurve(labels,avg(:,sq),1);
    Auc(sq) = auc; 
end
plot(sqenc_list,Auc,'-o')
ylim([0.5 0.85])


avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1),2));
NRT = 200;
labels = (1:NRT)>(NRT/2);
[~,~,~,auc] = perfcurve(labels,avg,1);
1-auc
