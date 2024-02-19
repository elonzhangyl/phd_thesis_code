clc
clear
%% Load data set:
load('und_seed101To320_1e6.mat')
dataset1 = dataset1(120*20+1:end,:,:);
dataset1 = cat(3,dataset1(1:end/2,:,:),dataset1(end/2+1:end,:,:));

% load('dam_scour0.24m_seed201To250_1e6.mat')
load('dam_scour0.36m_seed201To300_1e6.mat')
dataset2 = dataset2(120*20+1:end,:,:);
dataset21 = cat(3,dataset2(1:end/2,:,:),dataset2(end/2+1:end,:,:));

node = 1;
resam = 1;
perid = 9;
color_boundary = 0.3;

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
wd = 2048/4;
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
map = brewermap(3,'Set1');

figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.33]*36.36);
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
yticks(-4:2:2)
yticklabels({'10^{-4}','10^{-2}','10^{0}',...
    '10^{2}'})
xlabel('Frequency (Hz)')
ylabel('PSD ((m/s^2)^2/Hz)')
legend1 = legend([p1,p2],'Undamaged state','0.36 m damaged state',...
        'location','south','Orientation','vertical');
legend boxoff
set(legend1,...
    'Position',[0.436912467697463 0.227559866657652 0.31081450262083 0.156731827134704],...
    'FontSize',7);
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig3.PSD_boundary_owt.eps','Resolution',1000)



