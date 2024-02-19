clc
clear
%% Load data set:
load('und_seed101To112.mat')
load('dam_scour5m_seed201To225.mat')

% plot(dataset1(:,1,2))
% 
% wd = 256;
% [pxx1,f1] = pwelch(dataset1(:,1,2),wd,[],[],20);
% plot(f1,log10(pxx1))

% 
% plot(f1,log10(pxx1),f2,log10(pxx2),f3,log10(pxx3))
% legend('1','2','3')
% 
% spectrogram(dataset(:,5,93),[],[],[],320,'yaxis')
% 
% qqplot(dataset1(1:end/10,6,2))
%% sequence density ratio KL divergence for localization
dataset_ref = dataset1(1:7:end,6,1);
dataset_test = cat(3,dataset1(1:7:end,6,2:end),dataset2(1:7:end,6,1:11));

sqenc_list = 60;
alpha_list = 0;
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;%seed
NRR = size(dataset_ref,3);
NRT = size(dataset_test,3); 
DD = 1;
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
tic
for dd = 1:DD
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
                            X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                            X2 = zeros(size(Dataset200,1)-sqenc,sqenc);
                            for j=1:(size(Dataset100,1)-sqenc)
                                X1(j,:) = Dataset100(j:j+sqenc-1,ss);
                            end
                            for j=1:(size(Dataset200,1)-sqenc)
                                X2(j,:) = Dataset200(j:j+sqenc-1,ss);
                            end
                            x_de=X1';
                            x_nu=X2';
                            [wh_x_de,~] = KLIEP(x_de,x_nu);
                            SCORE_SS(ss) = 1/size(x_de,2)*sum(log(wh_x_de))
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
save('results_resample_60.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
    'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');

plot(SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:),'.');

% score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
% plot(score_avg(:),'.');


labels = (1:NRT)>11;

[~,~,~,auc] = perfcurve(labels,...
    abs(SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:)),1)







