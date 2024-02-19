clear;clc;


tic; 
fs = 20; ts = 1200; ns = ts*fs;%number of sample in one record
load('und_windwave_seed102To168_slt.mat');
dataset_ref = dataset1(1:ns*10,:);
% dataset_ref = reshape(dataset_ref,ns,10,4);
dataset_ref = reshape(dataset_ref,ts*20,10,4);
dataset_ref = dataset_ref(1:1:end,1:10,:);

load('dam_0.02_windwave_seed401To448_slt.mat');
dataset_test_2 = dataset1(1:ns*19,:);
% dataset_test_2 = reshape(dataset_test_2,ns,10,4);
dataset_test_2 = reshape(dataset_test_2,ts*20,19,4);
dataset_test_2 = dataset_test_2(1:1:end,1:19,:);

y = dataset_ref(:,1,4);parcorr(y,200);hold on
z = dataset_test_2(:,1,4);parcorr(z,200);

y = dataset_ref(:,1,4);autocorr(y,20);hold on
y = dataset_ref(:,2,4);autocorr(y,20);hold on
y = dataset_ref(:,3,4);autocorr(y,20);hold on
y = dataset_ref(:,4,4);autocorr(y,20);

y = dataset_ref(:,1,4);autocorr(y,20);hold on
z = dataset_test_2(:,1,4);autocorr(z,20);
