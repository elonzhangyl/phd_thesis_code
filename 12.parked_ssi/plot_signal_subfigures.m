clc
clear
node = 1;
fs = 20;
load('und_seed101To200_1e6.mat')
dataset1 = dataset1(120*20+1:11*60*fs,node,1);

load('dam_scour0.36m_seed201To250_1e6.mat')
dataset21 = dataset2(120*20+1:11*60*fs,node,1);


alpha = 0.0565;
dataset_test_1 = dataset1(:,1);
Dataset10 = dataset_test_1 + alpha*randn(size(dataset_test_1))...
                    .*(ones(size(dataset_test_1)).*rms(dataset_test_1));
dataset_test_2 = dataset21(:,1);
Dataset2 = dataset_test_2 + alpha*randn(size(dataset_test_2))...
                    .*(ones(size(dataset_test_2)).*rms(dataset_test_2));
%% plot
figurewidth = 9; %cm
ratio = 0.3;
f = figure('Position',[10 10 figurewidth figurewidth*ratio]*36.36);
plot(1:length(Dataset10),Dataset10)
xlim([0,9*60*fs])
xticks(0:3*60*fs:9*60*fs);
xticklabels({'0','180','360','540'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-15,15])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration_owt_a.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*ratio]*36.36);
plot(1:length(Dataset2),Dataset2)
xlim([0,9*60*fs])
xticks(0:3*60*fs:9*60*fs);
xticklabels({'0','180','360','540'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-15,15])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration_owt_b.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*ratio]*36.36);
plot(1:10*fs,Dataset10(1:10*fs))
xlim([0,10*fs])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-15,15])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration_owt_c.eps','Resolution',1000)

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*ratio]*36.36);
plot(1:10*fs,Dataset2(1:10*fs))
xlim([0,10*fs])
xticks(0:2*fs:10*fs);
xticklabels({'0','2','4','6','8','10'})
xlabel('Time (s)')
ylabel('Acc. (m/s^2)')
ylim([-15,15])
% title('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig2.acceleration_owt_d.eps','Resolution',1000)




