clc;clear;
%% Tuning parameter
ar_order = 30;
CoefSel = 2:15;
level_list = {'dam_no_water_6000s_soil7.84.mat'};
% {'dam_soil_lev7.6_10min_seeds.mat',...
%     'dam_soil_lev7_10min_seeds.mat',...
%     'dam_soil_lev6.4_6000s.mat'};
alpha_list = [0 0.10 0.20 0.30];
alpha = 0.3;
ReSmp = 1;%10:5:40
sensor_list=1:4;
n_rep = 10;%repeat times for getting errorbar
tic;    
SCORE_rep_errorbar_0 = zeros(length(sensor_list),40,...
    n_rep,length(alpha_list),length(level_list));

dl = 1;
%% Load data set no wave - 5min
load('und_10min_no water.mat');
Dataset10 = dataset1(1:end/2,:);
% Dataset10 = dataset1(end/2+1:end,[1 5 8]);
load('und_no_water_6000s.mat');
Dataset1 = zeros(6000,4,20);
for i=1:20
    Dataset1(:,:,i)=dataset1(size(dataset1,1)/20*(i-1)+1 ...
        :size(dataset1,1)/20*i,:);
end

filename_dam = string(level_list(dl));
load(filename_dam);
Dataset2 = zeros(6000,4,20);
for i=1:20
    Dataset2(:,:,i)=dataset2(size(dataset2,1)/20*(i-1)+1 ...
        :size(dataset2,1)/20*i,:);
end

Dataset20 = cat(3,Dataset1,Dataset2);


AR_Coef1 = zeros(size(Dataset20,3),length(CoefSel));
for nr = 1:size(Dataset20,3) %num of records
    dataset1 = Dataset10(1:ReSmp:end,:); %resample
    %noise
    dataset1 = dataset1 + alpha*randn(size(dataset1))...
                    .*(ones(size(dataset1)).*rms(dataset1));    
    dataset1 = zscore(dataset1); %normlization
    %% AR Coefficients
    m1 = ar(dataset1(:,1),ar_order,'burg');%sensor 1
    AR_Coef1(nr,:) = m1.A(CoefSel);
end
%% PCA
AR_Coef1 = AR_Coef1 - mean(AR_Coef1);
[PCA_Coef1,score,latent] = pca(AR_Coef1);







