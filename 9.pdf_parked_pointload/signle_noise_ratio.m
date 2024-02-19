clear
clc
sqenc = 1;
ReSmp = 1;
%% Load data set no wave - 5min
load('und_seed201To205.mat');
Dataset10 = dataset1(:,1,1);
dataset1 = Dataset10(1:ReSmp:end,:); %resample
noise = 0*randn(size(dataset1))...
                    .*(ones(size(dataset1)).*rms(dataset1));
                
r = snr(dataset1,noise) %dB

