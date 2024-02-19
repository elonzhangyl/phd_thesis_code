clear;clc;

tic;  
load('und_no_water_real_6000s_seed211To215.mat');
dataset_ref=dataset1;
dataset_ref = reshape(dataset_ref,1200*20,25,4);
dataset_ref = dataset_ref(:,:,1);
load('und_no_water_real_6000s_seed201To205.mat');
dataset_test_1 = dataset1;
dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);
dataset_test_1 = dataset_test_1(:,:,1);
dataset_test = cat(2,dataset_ref,dataset_test_1);
dataset_test = zscore(dataset_test);

%% Tuning parameter
sqenc_list = 205:210;%[1 5:5:300];
% sqenc_list = 1:10:200;
SQ = length(sqenc_list);
alpha = 0;
NCV = 5;
NR = 50;
NSML = round(NR/NCV);
NBIG = NR-NSML;
rng(1);my_indices = randperm(NR);
dataset_test = dataset_test + alpha*randn(size(dataset_test))...
    .*(ones(size(dataset_test)).*rms(dataset_test,1));
%% Damage Index
SCORE_NBIG_NSML_NCV_SQ = zeros(NBIG,NSML,NCV,SQ);
parfor sq = 1:SQ
    sq
    sqenc = sqenc_list(sq);
    SCORE_NBIG_NSML_NCV = zeros(NBIG,NSML,NCV);
    for ncv=1:NCV % NUM fold
        Dataset1 = dataset_test...
            (:,my_indices([1:10*(ncv-1),10*ncv+1:50]));
        Dataset2 = dataset_test(:,my_indices(10*(ncv-1)+1:10*ncv));
        SCORE_NBIG_NSML = zeros(NBIG,NSML);
        for nsml = 1:NSML
            Dataset200 = Dataset2(:,nsml);
            X2 = zeros(size(Dataset200,1)-sqenc,sqenc);
            for j=1:(size(Dataset200,1)-sqenc)
                X2(j,:) = Dataset200(j:j+sqenc-1);
            end
            SCORE_NBIG = zeros(NBIG,1);
            for nbig = 1:NBIG
                Dataset100 = Dataset1(:,nbig);
                X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
                for j=1:(size(Dataset100,1)-sqenc)
                    X1(j,:) = Dataset100(j:j+sqenc-1);
                end
                SCORE_NBIG(nbig) = KLDiv(X1,X2);  %KL(und,dam)
            end
            SCORE_NBIG_NSML(:,nsml) = SCORE_NBIG;
        end
        SCORE_NBIG_NSML_NCV(:,:,ncv) = SCORE_NBIG_NSML;
    end
    SCORE_NBIG_NSML_NCV_SQ(:,:,:,sq) = SCORE_NBIG_NSML_NCV;
end

mytime = toc;
save('sequence_sensitity_cv_0.15_test.mat','SCORE_NBIG_NSML_NCV_SQ',...
    'sqenc_list','SQ','alpha','NCV','NR','NSML','NBIG');



% SCORE_avg = squeeze(mean(SCORE_NBIG_NSML_NCV_SQ,1));
% SKL = zeros(NCV,SQ);
% for sq = 1:SQ
%     for ncv = 1:NCV
%         skl = SCORE_avg(:,ncv,sq);
% %         std(skl)
%         Skl = skl
%         SKL(ncv,sq) = sum(abs(skl - mean(skl)) > 2*std(skl))/NSML; 
%     end
% end
% aa=mean(SKL,1);
% plot(1:length(aa),aa)
% ylim([0 0.1]);

% SCORE_avg = squeeze(mean(SCORE_NBIG_NSML_NCV_SQ,1));
% SKL = zeros(SQ,1);
% for sq = 1:SQ
%     skl = SCORE_avg(:,:,sq);
% %         std(skl)
%     Skl = skl(:);
%     SKL(sq) = sum(abs(Skl - mean(Skl)) > 1.2*std(Skl))/NSML; 
% end
% aa=mean(SKL,1);
% plot(1:length(aa),aa)
% ylim([0 0.1]);

