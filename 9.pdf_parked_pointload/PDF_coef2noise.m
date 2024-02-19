clear;clc
ReSmp = 1;
nr = 22;
alpha = 0.30;
sqenc = 3;
sensor_list=1;

load('und_10min.mat');
Dataset10 = dataset1(1:end/2,[1 5 8]);
load('und_10min_seeds.mat');
Dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));
load('dam_soil_lev7.6_10min_seeds.mat');
Dataset20 = cat(3,Dataset1,Dataset2);


dataset1 = Dataset10(1:ReSmp:end,:); %resample
dataset2 = Dataset20(1:ReSmp:end,:,nr); %resample
%noise
dataset1 = dataset1 + alpha*randn(size(dataset1)).*dataset1;
dataset2 = dataset2 + alpha*randn(size(dataset2)).*dataset2;
dataset1 = zscore(dataset1); %normlization
dataset2 = zscore(dataset2); %normlization

X1=zeros(size(dataset1,1)-sqenc,sqenc);
X2=zeros(size(dataset2,1)-sqenc,sqenc);
for ss=1:length(sensor_list) % number of sensors
    sensor = sensor_list(ss);
    for j=1:(size(dataset1,1)-sqenc)
        X1(j,:) = (dataset1(j:j+sqenc-1,ss));
        X2(j,:) = (dataset2(j:j+sqenc-1,ss));
    end
    x_de=X1;
    x_nu=X2;
    pd_nu = mvnpdf(x_de,mean(x_nu),cov(x_nu));
    pd_de = mvnpdf(x_de,mean(x_de),cov(x_de));
    wh_x_de = 1/length(x_de)*sum(log(pd_nu./pd_de));% KL Div
%     SCORE_rep_errorbar_0(ss,nr,te,ns,dl) = wh_x_de;
end
wh_x_de