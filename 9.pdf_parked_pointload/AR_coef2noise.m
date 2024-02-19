clear;clc
load('und_10min.mat');
Dataset10 = dataset1(1:end/2,[1 5 8]);
load('und_10min_seeds.mat');
Dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));
load('dam_soil_lev7.6_10min_seeds.mat');
Dataset20 = cat(3,Dataset1,Dataset2);

% rng(1)
dataset1 = Dataset10; %resample
dataset1 = dataset1 + 0.3*randn(size(dataset1)).*dataset1;%noise
dataset1 = zscore(dataset1); %normlization


m1 = ar(dataset1(:,1),22,'burg');
AR_Coef1 = m1.A(2:22)'