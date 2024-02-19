clc;clear;
%% Load data set no wave - 5min
load('duffingCrack_dam1_amp1_seed1001_1200.mat')
dataset1 = dataset2;
load('duffingCrack_dam0.99_amp1_seed99001_99120.mat')
dataset21 = dataset2;

alpha = 0.032;
dataset_test_1 = dataset1(:,1);
Dataset10 = dataset_test_1 + alpha*randn(size(dataset_test_1))...
                    .*(ones(size(dataset_test_1)).*rms(dataset_test_1));
dataset_test_2 = dataset21(:,1);
Dataset2 = dataset_test_2 + alpha*randn(size(dataset_test_2))...
                    .*(ones(size(dataset_test_2)).*rms(dataset_test_2));
fs = 512;
%% plot
figurewidth = 19; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.33]*29.1);subplot(2,2,1)
plot(1:length(Dataset10),Dataset10)
xlim([0,5120])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
ylim([-1.3,1.3])
title('(a)')

subplot(2,2,2)
plot(1:length(Dataset2),Dataset2)
xlim([0,5120])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylim([-1.3,1.3])
title('(b)')
% ylabel('Acceleration signals at Sensor 1 (m/s^2)')

subplot(2,2,3)
plot(1:200,Dataset10(1:200,1))
xticks(0:40:200)
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
ylim([-1.3,1.3])
title('(c)')

subplot(2,2,4)
plot(1:200,Dataset2(1:200,1))
xticks(0:40:200)
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylim([-1.3,1.3])
title('(d)')

set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration.eps','Resolution',1000)

% get(figure,'default')

% set(findall(gcf,'-property','linewidth'),'linewidth',1)
% [acf,lags,~] = autocorr(Dataset10,'NumLags',20);
% plot(lags,acf,'-o')

%% for presentation
load('und_no_water_real_6000s_seed201To205.mat');
Dataset10 = dataset1(3:101,1);
figure('position',[100 100 800 120])
p = plot(1:99,Dataset10);
p.LineStyle = '-';
% p.LineWidth = 1;
p.Marker = 'o';
p.MarkerSize = 4;
p.MarkerEdgeColor = [0, 0.4470, 0.7410];
p.MarkerFaceColor = 'w';
xlabel('Time');
ylabel('Acceleration')
xticks([])
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% 2d pdf
[X,Y] = meshgrid(-5:0.1:5,-5:0.1:5);
Z = zeros(length(X),length(Y));
for i = 1:length(X)
    for j = 1:length(Y)
        Z(i,j) = mvnpdf([X(i,j),Y(i,j)],[0 0],[1 0.6; 0.6 2]) +...
            5*mvnpdf([X(i,j),Y(i,j)],[1 0.1],[1 0.6; 0.6 2]);
    end
end
surf(X,Y,Z)
grid off
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'ZTick',[]);

%1d pdf
x = -4:0.01:5;
y = normpdf(x,0,1)+normpdf(x,1,1)+2*normpdf(x,2,1);
plot(x,y);
grid off
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'ZTick',[]);

load('und_no_water_real_6000s_seed201To205.mat');
Dataset10 = dataset1(1:6000,1);
figure('position',[100 100 300 100])
plot(1:length(Dataset10),Dataset10)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'ZTick',[]);

Dataset10 = dataset1(6000+1:6000+6000,1);
figure('position',[100 100 300 100])
plot(1:length(Dataset10),Dataset10)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'ZTick',[]);




