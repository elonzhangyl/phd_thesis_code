clc
clear
%% Load data set:
load('und_seed101_1e3_chirp.mat')
x = dataset1(:,14);
% load('und_seed101To112.mat')
% x = dataset1(:,1,1);
% information
% get(0,'ScreenSize') and = [1 1 1536 864]
% screen physical size 20.75 inch * 11.67 inch
% screen: 74 pixels per inch; 29.1 pixels per centimeter
% 3 columns 6.3cm; 2 columns 9cm; 1 column 18cm; 1.5 column 14cm

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.78]*29.1);
spectrogram(x,256*2,50,[],20,'yaxis');
% ylim([0,1.5]);
hold on
% for i = [1 3 5 7 9 11 13]
%     xln = 0:0.01:20;
%     yln = (0.2+0.1/20*xln)*i;
%     plot(xln,yln,'--r','linewidth',0.5);
%     text(xln(end)-1,yln(end),[num2str(i),'x'],...
%         'HorizontalAlignment','right');
%     hold on
% end
hold off
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'spectrogram_chirp_1e3.eps','Resolution',1000)





