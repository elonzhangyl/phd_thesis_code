clc
clear
%% Load data set:
load('und_seed201To205.mat');
dataset_test_1 = dataset1;
dataset_test_1 = reshape(dataset_test_1,1200*20,25,4);
% load('dam_0.02_seed301To305.mat');
load('dam_0.04_seed401To405.mat');
dataset_test_2 = dataset2;
dataset_test_2 = reshape(dataset_test_2,1200*20,25,4);
color_boundary = 0.3;
map = brewermap(3,'Set1');
wd = 4096/4;
fs = 20;
nfft = 2^nextpow2(wd);
pwelch(dataset_test_1(:,1,:),wd,[],nfft,fs)

%% sequence density ratio KL divergence for localization
node = 4;
dataset_test = cat(3,dataset_test_1(:,1:25,node),dataset_test_2(:,1:25,node));

x = dataset_test;
wd = 4096/4;
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

figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.4]*36.36);
for i = 1
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p1 = plot(xplot, pxxlog_mean(:,i), 'color',map(2,:),...
    'LineWidth', 1);
hold on
p11 = plot(xplot, curve1, '-','color',map(2,:),...
    'LineWidth', 1);
p11.Color(4) = color_boundary;
hold on
p12 = plot(xplot, curve2, '-','color',map(2,:),...
    'LineWidth', 1);
p12.Color(4) = color_boundary;
end

for i = 2
xplot = f1';
curve1 = pxxlog_mean(:,i)'-3*pxxlog_std(:,i)';
curve2 = pxxlog_mean(:,i)'+3*pxxlog_std(:,i)';
p2 = plot(xplot, pxxlog_mean(:,i), '--','color',map(1,:),...
    'LineWidth', 1);
hold on
p21 = plot(xplot, curve1, '--','color',map(1,:),...
    'LineWidth', 0.5);
p21.Color(4) = color_boundary;
hold on
p22 = plot(xplot, curve2, '--','color',map(1,:),...
    'LineWidth', 0.5);
p22.Color(4) = color_boundary;
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
% yticklabels({'10^{-6}','10^{-4}','10^{-2}','10^{0}',...
%     '10^{2}'})
yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}',...
    '10^{0}','10^{1}','10^{2}','10^{3}'})
yticks([-4:3])
ylim([-4 3])
xlabel('Frequency (Hz)')
ylabel('PSD ((m/s^2)^2/Hz)')
title(strcat('Sensor',num2str(node)))
legend([p1,p2],'Undamaged','4% stiffness reduction',...
        'location','southeast','Orientation','vertical')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig3.PSD.eps','Resolution',1000)



