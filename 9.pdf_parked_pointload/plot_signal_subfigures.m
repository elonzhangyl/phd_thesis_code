clc;clear;
%% Load data set no wave - 5min
load('und_seed201To205.mat');
dataset11 = dataset1;
load('dam_0.02_seed301To305.mat');
dataset21 = dataset2;

alpha = 0.0315;
dataset_test_1 = dataset11(:,1);
Dataset10 = dataset_test_1 + alpha*randn(size(dataset_test_1))...
                    .*(ones(size(dataset_test_1)).*rms(dataset_test_1));
Dataset10 = zscore(Dataset10)/3;                
dataset_test_2 = dataset21(:,1);
Dataset2 = dataset_test_2 + alpha*randn(size(dataset_test_2))...
                    .*(ones(size(dataset_test_2)).*rms(dataset_test_2));
Dataset2 = zscore(Dataset2)/3;              
fs = 20;
%% plot
figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.3]*36.36);
plot(1:length(Dataset10),Dataset10)
xlim([0,12000])
xticks(0:2*60*fs:10*60*fs);
xticklabels({'0','120','240','360','480','600'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig.acceleration_a.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.3]*36.36);
plot(1:length(Dataset2),Dataset2)
xlim([0,12000])
xticks(0:2*60*fs:10*60*fs);
xticklabels({'0','120','240','360','480','600'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig.acceleration_b.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.3]*36.36);
plot(1:10*fs,Dataset10((1:10*fs)+fs*300,1))
xlim([0,10*fs])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
% title('(c)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig.acceleration_c.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.3]*36.36);
plot(1:10*fs,Dataset2((1:10*fs)+fs*300,1))
xticks(0:40:200)
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-1.4,1.4])
% title('(d)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig.acceleration_d.eps','Resolution',1000)

% get(figure,'default')

% set(findall(gcf,'-property','linewidth'),'linewidth',1)
% [acf,lags,~] = autocorr(Dataset10,'NumLags',20);
% plot(lags,acf,'-o')



