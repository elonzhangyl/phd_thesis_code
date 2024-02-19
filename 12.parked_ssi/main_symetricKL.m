clc
clear
%% Load data set:
load('und_seed101To150.mat')
% load('dam_scour0.6m_seed201To210.mat')
% load('dam_scour1.2m_seed201To210.mat')
% load('dam_scour1.8m_seed201To210.mat')
% load('dam_scour2.4m_seed201To210.mat')
% load('dam_scour3.6m_seed201To210.mat')
load('dam_scour4.8m_seed201To210.mat')
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
dataset_ref = dataset1(:,5,21);
dataset_test = cat(3,dataset1(:,5,11:20),dataset2(:,5,1:10));

wd_list = 1;%[0.7,0.8,0.9,1];
sqenc_list = 20;%5:5:80;
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
                    Dataset2 = dataset_test(:,:,nrt); %choose record 
    %                 rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = (Dataset20); %normlization 
                    SCORE_SS_NRR = zeros(SS,NRR);
                    for nrr = 1:NRR %NUM RECORD reference
                        Dataset1 = dataset_ref(:,:,nrr); %choose record
                         %choose record
    %                     rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = (Dataset10); %normlization
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
                            if nrt == 1
%                                 [~,wh_x_re,x_ce,alphah] = KLIEP2(x_de,x_nu,x_nu,wd);
%                                 SCORE_SS(ss) = mean(log(wh_x_re))
%                                SCORE_SS(ss) = KLDiv(X1,X2)
                                [~,wh_x_re1,x_ce1,alphah1] = KLIEP2(x_de,x_nu,x_nu,wd); 
                                [~,wh_x_re2,x_ce2,alphah2] = KLIEP2(x_nu,x_de,x_de,wd);              

                                SCORE_SS(ss) = mean(log(wh_x_re1))+mean(log(wh_x_re2))

                            else
                                X_nu1=kernel_Gaussian(x_nu,x_ce1,wd);
                                wh_x_re1=(X_nu1*alphah1)';
                                X_de2=kernel_Gaussian(x_de,x_ce2,wd);
                                wh_x_re2=(X_de2*alphah2)';
                                SCORE_SS(ss) = mean(log(wh_x_re1))+mean(log(wh_x_re2))
                            end

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
plot(SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:),'.');
% 
% % score_avg = mean(SCORE_SS_NRR_NRT_SN_NL_SQ_DD,2);
% % plot(score_avg(:),'.');
% 
% 
% labels = (1:NRT)>10;
% 
% [~,~,~,auc] = perfcurve(labels,...
%     (SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:)),1)



labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NRR,SN,NL,SQ,DD);

for dd = 1:DD
    for sq = 1:SQ
        for nl = 1:NL % noise level
            for sn = 1:SN % seed noise
                for nrt=1:NRT % NUM RECORD TEST
                    for nrr = 1:NRR %NUM RECORD reference
                        for ss = 1:SS
                            [~,~,~,auc] = perfcurve(labels,...
                                squeeze(SCORE_SS_NRR_NRT_SN_NL_SQ_DD...
                                (ss,nrr,:,sn,nl,sq,dd)),1);
                            AUC(ss,nrr,sn,nl,sq,dd) = auc;
                        end
                    end
                end
            end
        end
    end
end 
squeeze(AUC)



