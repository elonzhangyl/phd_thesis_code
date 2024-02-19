clc
clear
%% Load data set:
load('und_seed101To320_1e6.mat')
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

load('dam_scour0.6m_seed201To260_1e6.mat')
dataset2 = dataset2(120*20+1:end,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));

node = 1;
resam = 1;
perid = 9;

dd=1;
alpha = 0.032;
dataset_test_1 = squeeze(dataset1(1:resam:20*perid*60,node,[61:70,181:190]));
dataset_test_1 = dataset_test_1 + alpha*randn(size(dataset_test_1))...
                    .*(ones(size(dataset_test_1)).*rms(dataset_test_1));
dataset_test_2 = squeeze(dataset21(1:resam:20*perid*60,node,[1:10,61:70]));
dataset_test_2 = dataset_test_2 + alpha*randn(size(dataset_test_2))...
                    .*(ones(size(dataset_test_2)).*rms(dataset_test_2));
%% sequence density ratio KL divergence for localization
dataset_test = cat(3,dataset_test_1,dataset_test_2);

x = dataset_test;
wd = 2048;
fs = 20;
nfft = 2^nextpow2(wd);
%  [pxx,f1] = pwelch(x(:,1,1),wd,[],[],fs);
pxx = zeros( (nfft/2 + 1) ,size(x,2),size(x,3));
for i = 1:size(x,2)
    for j = 1:size(x,3)
        [pxx(:,i,j),f1] = pwelch(x(:,i,j),wd,[],nfft,fs);
    end
end
pxxlog = log10(pxx);

pxxlog_mean = squeeze(mean(pxxlog,2));
pxxlog_std = squeeze(std(pxxlog,[],2));
map = brewermap(3,'Set1');
figurewidth = 19; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.33]*29.1);
for i = 1
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p1 = plot(xplot, pxxlog_mean(:,i), 'color',map(2,:),...
    'LineWidth', 1,'HandleVisibility','off');
hold on
x2 = [xplot, fliplr(xplot)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween,map(2,:),'FaceAlpha',0.5,'linestyle','none','edgecolor','none');
end

for i = 2
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p2 = plot(xplot, pxxlog_mean(:,i), '.','color',map(1,:),...
    'markersize', 3,'HandleVisibility','off');
hold on
x2 = [xplot, fliplr(xplot)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, map(1,:),'FaceAlpha', 0.5,'linestyle','none','edgecolor','none');
end

% for i = 3
% xplot = f1';
% curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
% curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
% plot(xplot, pxxlog_mean(:,i), 'g', 'LineWidth', 1);
% hold on
% x2 = [xplot, fliplr(xplot)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g','FaceAlpha', 0.5,'linestyle','none');
% end

% xticklabels({'0','2','4','6','8','10'})
yticklabels({'10^{-6}','10^{-4}','10^{-2}','10^{0}',...
    '10^{2}'})
xlabel('Frequency (Hz)')
ylabel('PSD ((m/s^2)^2/Hz)')
legend('Undamaged','0.60 m damaged',...
        'location','northeast','Orientation','vertical')
legend boxoff
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig3.psd_shaded.eps','Resolution',1000)



