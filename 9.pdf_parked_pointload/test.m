clear;clc;

tic;
%% Tuning parameter
sqenc_list = 150;
alpha_list = 0.15;
sensor_list = 1;
SQ = length(sqenc_list);
NL = length(alpha_list);
SS = length(sensor_list);
SN = 1;
NRR = 25;
NRT = 50; 
DD = 1;
SCORE_SS_NRR_NRT_SN_NL_SQ_DD = zeros(SS,NRR,NRT,SN,NL,SQ,DD);
for dd = 1:DD
    dd
    file_ref = {'und_seed201To205.mat',...
        'und_seed211To215.mat',...
        'und_seed221To225.mat',...
        'und_seed231To235.mat',...
        'und_seed241To245.mat',...
        'und_seed201To205.mat',...
        'und_seed211To215.mat',...
        'und_seed221To225.mat',...
        'und_seed231To235.mat',...
        'und_seed241To245.mat'};
    load(string(file_ref(dd)));
    dataset_ref=dataset1;
    dataset_ref = reshape(dataset_ref,1200*20,25,4);
    
    file_test_1 = {'und_seed206To210.mat',...
        'und_seed216To220.mat',...
        'und_seed226To230.mat',...
        'und_seed236To240.mat',...
        'und_seed246To250.mat',...
        'und_seed206To210.mat',...
        'und_seed216To220.mat',...
        'und_seed226To230.mat',...
        'und_seed236To240.mat',...
        'und_seed246To250.mat'};
    load(string(file_test_1(dd)));
    dataset_test_1 = dataset1;
    dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);


%     file_test_2 = {'dam_0.02_seed301To305.mat',...
%         'dam_0.02_seed306To310.mat',...
%         'dam_0.02_seed311To315.mat',...
%         'dam_0.02_seed316To320.mat',...
%         'dam_0.02_seed321To325.mat',...
%         'dam_0.04_seed401To405.mat',...
%         'dam_0.04_seed406To410.mat',...
%         'dam_0.04_seed411To415.mat',...
%         'dam_0.04_seed416To420.mat',...
%         'dam_0.04_seed421To425.mat'};
%     load(string(file_test_2(dd)));
    
    file_test_2 = {'dam_0.01_seed1001To1005.mat',...
        'dam_0.01_seed1006To1010.mat',...
        'dam_0.01_seed1011To1015.mat',...
        'dam_0.01_seed1016To1020.mat',...
        'dam_0.01_seed1021To1025.mat',...
        'dam_0.03_seed3001To3005.mat',...
        'dam_0.03_seed3006To3010.mat',...
        'dam_0.03_seed3011To3015.mat',...
        'dam_0.03_seed3016To3020.mat',...
        'dam_0.03_seed3021To3025.mat'};
    load(string(file_test_2(dd)));

    dataset_test_2 = dataset2;
    dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);
    dataset_test = cat(2,dataset_test_1,dataset_test_2);

  WW=1;
    %% Damage Index
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
                            SCORE_WW(ww) = KLDiv(X1,X2);%de in front
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
% save('result_pdf_KLDiv_k155.mat','SCORE_SS_NRR_NRT_SN_NL_SQ_DD',...
%     'SS','NRT','NL','SN','SQ','DD','alpha_list','sqenc_list','sensor_list');
asss = squeeze(SCORE_WW_NRR_NRT_SN_NL_SQ_DD);


