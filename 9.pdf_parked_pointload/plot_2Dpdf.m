clear
clc
sqenc = 50;
load('und_no_water_real_6000s_seed201To205.mat');
dataset1 = dataset1 + 0*randn(size(dataset1))...
            .*(ones(size(dataset1)).*rms(dataset1));
dataset1 = zscore(dataset1); %normlization
X1=zeros(size(dataset1,1)-sqenc,sqenc);
for ss=1 %sensors
    for j=1:(size(dataset1,1)-sqenc)
        X1(j,:) = (dataset1(j:j+sqenc-1,ss));
    end
end
% cov(X1)
X1 = X1(1:12001,:);
[acf1,~,~] = autocorr(X1(:,1),100);


load('dam_no_water_real_allmatrix_soil7.20_seed301.mat');
dataset2 = dataset2 + 0*randn(size(dataset2))...
            .*(ones(size(dataset2)).*rms(dataset2));
dataset2 = zscore(dataset2); %normlization
X2=zeros(size(dataset2,1)-sqenc,sqenc);
for ss=1 %sensors
    for j=1:(size(dataset2,1)-sqenc)
        X2(j,:) = (dataset2(j:j+sqenc-1,ss));
    end
end
% cov(X2)
[acf2,~,~] = autocorr(X2(:,1),100);
plot(1:101,[acf1-acf2])
scatter(X1(:,1),X1(:,44),'.');hold on
scatter(X2(:,1),X2(:,44),'.')



% scatter(X1(:,1),X1(:,2),'.')
% hist3(X1,'CDataMode','auto')
nbins=[20 20];
for i = 1:sqenc
    subplot(1,3,i)
    if i==1 || i==2
    [N,C]=hist3(X1(:,i:i+1),nbins);
    end
    if i==3
        [N,C]=hist3(X1(:,[1,3]),nbins);
    end
    contourf(C{1},C{2},N)
    axis equal
    xlim([-3 3])
    ylim([-3 3])
    xlabel('X_1')
end
%     xlabel('$x_k$','interpreter','latex')
%     ylabel('$x_{k+1}$','interpreter','latex')
for i = 1:sqenc
    subplot(5,6,i)
    [N,C]=hist3(X1(:,[1,i]),nbins);
    contourf(C{1},C{2},N)
    axis equal
    xlim([-3 3])
    ylim([-3 3])
    xlabel('X_1')
end

