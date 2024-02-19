clc
clear
%% Load data set:
load('und_seed101To112.mat')
% load('dam_scour0.12m_seed201To207.mat')
load('dam_scour0.12m_seed201To210.mat')
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
node = 14;
% dataset11 = zeros(size(dataset1));
% for i = 1:10
% dataset11(:,node,i) = dataset1(:,node,i) + 0.3*randn(size(dataset1(:,node,i)))...
%                         .*(ones(size(dataset1(:,node,i))).*rms(dataset1(:,node,i),1));
% std(dataset11(:,node,i));
% end
% dataset1 = dataset11;

dataset_ref = dataset1(:,:,1);
% dataset_test = cat(3,dataset1(:,node,11:20),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10),...
%     dataset25(:,node,1:10),dataset26(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10),dataset24(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,2:11),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10),dataset23(:,node,1:10));
% dataset_test = cat(3,dataset1(:,node,1:10),dataset21(:,node,1:10),...
%     dataset22(:,node,1:10));
dataset_test = cat(3,dataset1(:,:,2:11),dataset21(:,:,1:10));

% aa = std(squeeze(dataset_test))';
% plot(aa,'.')
% xx = [dataset_test(:,:,1),dataset_test(:,:,31),dataset_test(:,:,41)];
% plot(xx(:))
% std(dataset_test(:,:,11))
% 
% for i = 1:size(dataset_test,3)
%     a(i) = std(dataset_test(:,1,i));
% end
% plot(a,'.')
% 
% load('dam_scour0.12m_seed203To203.mat')
% d1 = dataset2(:,14,1);
% load('dam_scour0.24m_seed203To203.mat')
% d2 = dataset2(:,14,1);
% load('dam_scour0.36m_seed203To203.mat')
% d3 = dataset2(:,14,1);
% std(d1)
% std(d2)
% std(d3)

% a(11:20) = a(11:20)-(1.745-std(d1));plot(a,'.')
%% delete effect of wind speed 0.1m/s
% dataset_ref =dataset_ref(100*20+1:end,:,:);
% dataset_test =dataset_test(100*20+1:end,:,:);

%% main
wd_list = 1;%[0.7,0.8,0.9,1];
alpha_list = 0;
NL = length(alpha_list);
ss = [12 1 5 8 11];
NRR = size(dataset_ref,3);
NRT = size(dataset_test,3); 
DD = length(wd_list);
SCORE_NRR_NRT_NL_DD = zeros(NRR,NRT,NL,DD);
tic
% mu = zeros(sqenc_list,1);
% Sigma = eye(sqenc_list) * cov(dataset_ref)*0.36;
% x_re = mvnrnd(mu,Sigma,20000)';
for dd = 1:DD
    wd = wd_list(dd);
    SCORE_NRR_NRT_NL = zeros(NRR,NRT,NL);
        for nl = 1:NL % noise level
%             nla
            alpha =alpha_list(nl);
            SCORE_NRR_NRT = zeros(NRR,NRT);
            for nrt=1:NRT % NUM RECORD TEST
                nrt
                Dataset2 = dataset_test(:,:,nrt); %choose record 
                Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                    .*(ones(size(Dataset2)).*rms(Dataset2,1));
                Dataset200 = zscore(Dataset20); %normlization 
                SCORE_NRR = zeros(NRR,1);
                for nrr = 1:NRR %NUM RECORD reference
                    Dataset1 = dataset_ref(:,:,nrr); %choose record
                    Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                            .*(ones(size(Dataset1)).*rms(Dataset1,1));
                    Dataset100 = zscore(Dataset10); %normlization
                    X1 = Dataset100(:,ss);
                    X2 = Dataset200(:,ss);
                    x_de=X1';
                    x_nu=X2';
                    [~,wh_x_re,~,~] = KLIEP2(x_de,x_nu,x_nu,wd);  
                    SCORE_NRR(nrr) = mean(log(wh_x_re))
%                                SCORE_SS(ss) = KLDiv(X1,X2)
%                             SCORE_SS(ss) = std(X2(:,1))
                end
                SCORE_NRR_NRT(:,nrt) = SCORE_NRR; 
            end
        SCORE_NRR_NRT_NL(:,:,nl) = SCORE_NRR_NRT;
        end
    SCORE_NRR_NRT_NL_DD(:,:,:,:,dd) = SCORE_NRR_NRT_NL;
end
mytime = toc;
% save('results_10m_sq80.mat','SCORE_NRR_NRT_NL_DD',...
%     'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');
% 
figure('position',[100 100 600 150])
plot(SCORE_NRR_NRT_NL_DD(:),'.','MarkerSize',20);
% % % 
% % % score_avg = mean(SCORE_NRR_NRT_NL_DD,2);
% % % plot(score_avg(:),'.');
% % 
% % 
% % labels = (1:NRT)>10;
% % 
% % [~,~,~,auc] = perfcurve(labels,...
% %     (SCORE_NRR_NRT_NL_DD(:)),1)
% 
% 
NRT = 20
labels = (1:NRT)>(NRT/2);
AUC = zeros(NRR,NL,DD);
for dd = 1:DD
    for nl = 1:NL % noise level
        for nrr = 1:NRR %NUM RECORD reference
                [~,~,~,auc] = perfcurve(labels,...
                    squeeze(SCORE_NRR_NRT_NL_DD...
                    (nrr,[1:10,11:20],nl,dd)),1);
                AUC(nrr,nl,dd) = auc;
        end
    end
end 
squeeze(AUC)



