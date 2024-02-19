clear
fs = 200;
t = 0:1/fs:10;
x = randn(length(t),1);

plot(t,x)
plot(t,x.^2-x)

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