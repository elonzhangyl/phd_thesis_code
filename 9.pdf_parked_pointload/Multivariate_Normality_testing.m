clear;clc
sqenc = 4;
ReSmp = 1;
%% Load data set no wave - 5min
% load('und_10min.mat');
% Dataset10 = dataset1(1:end/2,[1 5 8]);
% load('und_10min_no_water_real.mat');Dataset10= dataset1;

load('und_no_water_6000s.mat');Dataset10= dataset1(1:6000,:);

dataset1 = Dataset10(1:ReSmp:end,:); %resample
% dataset1 = dataset1 + 0.05*randn(size(dataset1)).*dataset1;% add noise
dataset1 = dataset1 + 0.0*randn(size(dataset1))...
                    .*(ones(size(dataset1)).*rms(dataset1));
dataset1 = zscore(dataset1); %normlization
X1=zeros(size(dataset1,1)-sqenc+1,sqenc);
for ss=2 %sensors
    for j=1:(size(dataset1,1)-sqenc+1)
        X1(j,:) = (dataset1(j:j+sqenc-1,ss));
    end
end
                    
[H stats] = mardiatest(X1,0.05)
% lillietest(X1(:,1))
% qqplot(X1(:,2))
% skewness(X1) 
% kurtosis(X1)
% xx = randn(120000,1);
% qqplot(xx)

% 
% 
% [N,c] = hist3([dataset1(1:120000-30,1,1),dataset1(1+30:120000,1,1)],...
%     'Nbins',[50 50]);
% c1 = c(:,1);c2 = c(:,2);
% data = [c1{:};c2{:}]';
% x=linspace(min(data(:,1)),max(data(:,1)),50);
% y=linspace(min(data(:,2)),max(data(:,2)),50);
% [X,Y]=meshgrid(x,y);
% F=interp2(data(:,1),data(:,2),N,X,Y);
% contourf(X,Y,F,'LineColor','none')
% 
% dataset1 = [dataset1(:,1,1);dataset1(:,1,2);dataset1(:,1,3);dataset1(:,1,4)...
%     ;dataset1(:,1,5);dataset1(:,1,6);dataset1(:,1,7);dataset1(:,1,8)...
%     ;dataset1(:,1,9);dataset1(:,1,10)];