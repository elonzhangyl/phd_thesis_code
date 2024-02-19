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
pca_order = 13;
sqenc_list = 130;%5:5:80;
alpha_list = 0.032;
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
for dd = 1
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
    for sq = 1
        sqenc = sqenc_list(sq);
        SCORE_WW_NRR_NRT_SN_NL = zeros(WW,NRR,NRT,SN,NL);
        for nl = 1 % noise level
%             nla
            alpha =alpha_list(nl);
            SCORE_WW_NRR_NRT_SN = zeros(WW,NRR,NRT,SN);
            for sn = 1 % seed noise
                SCORE_WW_NRR_NRT = zeros(WW,NRR,NRT);
                for nrt=1 % NUM RECORD TEST
                    nrt
                    Dataset2 = dataset_test(:,nrt); %choose record 
                    rng(100*(sn-1)+NRR+nrt);
                    Dataset20 = Dataset2 + alpha*randn(size(Dataset2))...
                        .*(ones(size(Dataset2)).*rms(Dataset2,1));
                    Dataset200 = (Dataset20); %normlization 
                    SCORE_WW_NRR = zeros(WW,NRR);
                    for nrr = 1 %NUM RECORD reference
                        Dataset1 = dataset_ref(:,nrr); %choose record
                         %choose record
                        rng(100*(sn-1)+nrr);
                        Dataset10 = Dataset1 + alpha*randn(size(Dataset1))...
                                .*(ones(size(Dataset1)).*rms(Dataset1,1));
                        Dataset100 = (Dataset10); %normlization
                        SCORE_WW = zeros(WW,1);

                        X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                        for j=1:(size(Dataset100,1)-sqenc)
                            X1(j,:) = Dataset100(j:j+sqenc-1);
                        end
                        [coeff,score,latent,~,explain] = pca(X1);
                        X1_PCA = score(:,1:pca_order);
                        X1_PCA = zscore(X1_PCA);
                        for ww = 1:WW % width
                            wd = wd_list(ww);
                            %% KLEIP PCA
                            x_de=X1_PCA';
                            distan = zeros(size(x_de,2),size(x_de,2));
                            for aa = 1:size(x_de,2)
                                aa
                                for bb = 1:size(x_de,2)
                                    distan(aa,bb) = norm(x_de(aa)-x_de(bb));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


mytime = toc;
median(distan(:))
mean(distan(:))
median(min(distan(distan>0)))
mean(min(distan(distan>0)))
median(max(distan))
mean(max(distan))
