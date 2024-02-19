clc;clear;
%% Load data set no wave - 5min
load('und_seed201To205.mat');
dataset11 = dataset1;
alpha = 0;
dataset_test_1 = dataset11(:,1);
Dataset10 = dataset_test_1 + alpha*randn(size(dataset_test_1))...
                    .*(ones(size(dataset_test_1)).*rms(dataset_test_1));
Dataset10 = zscore(Dataset10)/3;  

load('dam_0.04_seed401To405.mat');
dataset21 = dataset2;
dataset_test_2 = dataset21(:,1);
Dataset2 = dataset_test_2 + alpha*randn(size(dataset_test_2))...
                    .*(ones(size(dataset_test_2)).*rms(dataset_test_2));
Dataset2 = zscore(Dataset2)/3;              
fs = 20;  
%% plot
figurewidth = 19; %cm
f = figure('outerPosition',[10 10 figurewidth figurewidth*0.45]*36.36);
subplot(2,2,1)
plot(1:length(Dataset10),Dataset10)
xlim([0,24000])
xticks(0:200*fs:20*60*fs);
xticklabels({'0','200','400','600','800','1000','1200'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
title('(a)')

subplot(2,2,2)
plot(1:length(Dataset2),Dataset2)
xlim([0,24000])
xticks(0:200*fs:20*60*fs);
xticklabels({'0','200','400','600','800','1000','1200'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
title('(b)')
% ylabel('Acceleration signals at Sensor 1 (m/s^2)')

subplot(2,2,3)
plot(1:10*fs,Dataset10((1:10*fs)+fs*300,1))
xlim([0,10*fs])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
title('(c)')

subplot(2,2,4)
plot(1:10*fs,Dataset2((1:10*fs)+fs*300,1))
xticks(0:40:200)
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
title('(d)')

set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration.eps','Resolution',1000)


%% for presentation
% load('und_no_water_real_6000s_seed201To205.mat');
% Dataset10 = dataset1(3:101,1);
% figure('position',[100 100 800 120])
% p = plot(1:99,Dataset10);
% p.LineStyle = '-';
% % p.LineWidth = 1;
% p.Marker = 'o';
% p.MarkerSize = 4;
% p.MarkerEdgeColor = [0, 0.4470, 0.7410];
% p.MarkerFaceColor = 'w';
% xlabel('Time');
% ylabel('Acceleration')
% xticks([])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% % 2d pdf
% [X,Y] = meshgrid(-5:0.1:5,-5:0.1:5);
% Z = zeros(length(X),length(Y));
% for i = 1:length(X)
%     for j = 1:length(Y)
%         Z(i,j) = mvnpdf([X(i,j),Y(i,j)],[0 0],[1 0.6; 0.6 2]) +...
%             5*mvnpdf([X(i,j),Y(i,j)],[1 0.1],[1 0.6; 0.6 2]);
%     end
% end
% surf(X,Y,Z)
% grid off
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'ZTickLabel',[]);
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% set(gca,'ZTick',[]);
% 
% %1d pdf
% x = -4:0.01:5;
% y = normpdf(x,0,1)+normpdf(x,1,1)+2*normpdf(x,2,1);
% plot(x,y);
% grid off
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'ZTickLabel',[]);
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% set(gca,'ZTick',[]);
% 
% load('und_no_water_real_6000s_seed201To205.mat');
% Dataset10 = dataset1(1:6000,1);
% figure('position',[100 100 300 100])
% plot(1:length(Dataset10),Dataset10)
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'ZTickLabel',[]);
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% set(gca,'ZTick',[]);
% 
% Dataset10 = dataset1(6000+1:6000+6000,1);
% figure('position',[100 100 300 100])
% plot(1:length(Dataset10),Dataset10)
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'ZTickLabel',[]);
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% set(gca,'ZTick',[]);




