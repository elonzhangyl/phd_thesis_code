clear
clc
node = 14;

load('gwn_seed101_1e4.mat')
noise1 = noise(201:end);
load('gwn_seed101_1e6.mat')
noise2 = noise(201:end);
load('und_seed101To101_1e4.mat')
dataset11 = dataset1(:,node);
load('und_seed101To101_1e6.mat')
dataset12 = dataset1(:,node);

% plot(dataset11,'.')
% qqplot(dataset11)
% 
% qqplot(dataset12)



% load('Copy_of_gwn_seed101_1e0.mat')
% noise10 = noise(201:end);
% load('Copy_of_gwn_seed101_1e4.mat')
% noise20 = noise(201:end);
% load('und_seed101To101_linear.mat')
% dataset110 = dataset1(:,node);
% load('und_seed101To110.mat')
% dataset120 = dataset1(:,node,1);

% [pxx1,f1] = pwelch(x,wd,[],nfft,20);
% [pxx2,f2] = pwelch(noise10,wd,[],nfft,20);
% plot(f1,log10(pxx1./pxx2));
% 
% [pxx1,f1] = pwelch(x,wd,[],nfft,20);
% [pxx2,f2] = pwelch(noise1,wd,[],nfft,20);
% plot(f1,log10(pxx1./pxx2));hold on
% 
% % [pxx1,f1] = pwelch(dataset120,wd,[],nfft,20);
% % [pxx2,f2] = pwelch(noise20,wd,[],nfft,20);
% % plot(f1,log10(pxx1));
% 
% [pxx1,f1] = pwelch(dataset12,wd,[],nfft,20);
% [pxx2,f2] = pwelch(noise2,wd,[],nfft,20);
% plot(f1,log10(pxx1));


%% Nonlinearities Effects
figurewidth = 18; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.33]*29.1);

wd = 2048;
nfft = 2000;
[pxx1,f1] = pwelch(x,wd,[],nfft,20);
[pxx2,f2] = pwelch(noise1,wd,[],nfft,20);
plot(f1,log10(pxx1./pxx2),'linewidth',1);hold on
% xlim([20,140])

[pxx3,f3] = pwelch(dataset12,wd,[],nfft,20);
[pxx4,f4] = pwelch(noise2,wd,[],nfft,20);
plot(f3,log10(pxx3./pxx4),'linewidth',1)
% xlim([20,140])
legend('Excitation magnitude 1E0','Excitation magnitude 1E6',...
        'location','southeast','Orientation','vertical')
% yticklabels({'10^{-14}','10^{-12}','10^{-10}','10^{-8}',...
%    })
xlabel('Frequency (Hz)')
ylabel('FRF magnitude ((m/s^2)^2/N)')
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'frf.eps','Resolution',1000)


