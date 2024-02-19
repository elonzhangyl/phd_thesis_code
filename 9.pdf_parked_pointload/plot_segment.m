load('dam_0.02_windwave_seed401To448_slt.mat')

plot(1:1200*200*2,dataset1(1:1200*200*2,1));
for i = 1:2
    xline(1200*200*i);
    text(1200*200/6+1200*200*(i-1),-0.7,['seed10',num2str(i)])
end
xticks([1200*200*0,1200*200*1,1200*200*2])
xticklabels({'0','20','40'})
xlabel('Time(min)');
ylabel('Support Structure x-acceleration (m/s^2)')
title('Tower base acceleration time-series under different seeds')