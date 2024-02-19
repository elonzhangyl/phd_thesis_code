clc
clear
%% Load data set:
% load('und_seed101To150_1e6.mat')
load('und_seed101To320_1e6.mat')
sqenc = 1;
ReSmp = 1;
%% Load data set no wave - 5min
Dataset10 = dataset1(60*2*20:11*60*20,1,1);
% noise = 0.032*randn(size(Dataset10))...
%                     .*(ones(size(Dataset10)).*rms(Dataset10));
% noise = 0.0565*randn(size(Dataset10))...
%                     .*(ones(size(Dataset10)).*rms(Dataset10));
noise = 0.08*randn(size(Dataset10))...
                    .*(ones(size(Dataset10)).*rms(Dataset10));
                
r = snr(Dataset10+noise,noise) %dB

