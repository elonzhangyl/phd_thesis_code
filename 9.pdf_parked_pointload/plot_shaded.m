clc
clear
%% Load data set:
load('und_seed201To205.mat');
dataset_test_1 = dataset1;
dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);

load('dam_0.02_seed301To305.mat');
dataset_test_2 = dataset2;
dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);

%% sequence density ratio KL divergence for localization
node = 1;
dataset_test = cat(3,dataset_test_1(:,1:25,node),dataset_test_2(:,1:25,node));

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

% plot(f1,pxxlog_mean(:,1),'.',f1,pxxlog_mean(:,2),'.',f1,pxxlog_mean(:,3),'.')
% pxxtop = zeros(size(x,1),size(x,2));
% pxxbot = zeros(size(x,1),size(x,2));
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         pxxtop(i,j) = quantile(pxxlog(i,j,:),0.995);
%         pxxbot(i,j) = quantile(pxxlog(i,j,:),0.005);
%     end
% end

pxxlog_std = squeeze(std(pxxlog,[],2));

figurewidth = 19; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.33]*29.1);
for i = 1
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p1 = plot(xplot, pxxlog_mean(:,i), 'color',[0.3010 0.7450 0.9330],...
    'LineWidth', 2);
hold on
x2 = [xplot, fliplr(xplot)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.3010 0.7450 0.9330]	,'FaceAlpha', 1,'linestyle','none');
end

for i = 2
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p2 = plot(xplot, pxxlog_mean(:,i), '.','color',[0.8500 0.3250 0.0980],...
    'markersize', 6);
hold on
x2 = [xplot, fliplr(xplot)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.8500 0.3250 0.0980]	,'FaceAlpha', 0.5,'linestyle','none');
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
legend([p1,p2],'Undamaged','2% damaged',...
        'location','southeast','Orientation','vertical')
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'psd_shaded.eps','Resolution',1000)



