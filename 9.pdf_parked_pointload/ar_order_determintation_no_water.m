%for decide AR order
clc;clear;

load('und_seed201To205.mat');

Z = dataset1(1:20*60*20,1);
% Z = autocorr(Z,'NumLags',300-1);
% figure;
% subplot(2,1,1);autocorr(Z,'NumLags',size(Z,1)-1);
% subplot(2,1,2);parcorr(Z,'NumLags',100-1);
%estimate ARMA coefficients
p = 300;q = 1;

LOGL = zeros(p,q); % Initialize
PQ = zeros(p,q);
tic
parfor i = 1:p
    i
    mod = arima(i,0,0);
    [fit,~,logL] = estimate(mod,Z,'Display','off');
    LOGL(i,1) = logL;
    PQ(i,1) = i+1;
end
mytime = toc;
%Step 4: Calculate the BIC.
numObs = size(Z,1);
LOGL = reshape(LOGL,p*q,1);
PQ = reshape(PQ,p*q,1);
[aic,~] = aicbic(LOGL,PQ+1,numObs);
aic = reshape(aic,p,q);

save('result_ar_order.mat','aic');

%% PLOT
xx=1:length(aic);yy=aic;
figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.62]*29.1);
p = plot(xx,yy);
% p.Color = 'k';
p.LineStyle = '-';
p.LineWidth = 1;
p.Marker = 'o';
p.MarkerSize = 2;
% p.MarkerEdgeColor = 'k';
% p.MarkerFaceColor = 'k';
xlim([1,50]);
% ylim([0 0.12]);
% xticks([1,50:50:300]);
xlabel('Model order');
ylabel('AIC');
% ytickformat('percentage')
% yline(0.05,'r');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'figA1.AIC.eps','Resolution',1000)
%% validation
initial_arima_model = arima(191,1,1);
arima_model = estimate(initial_arima_model,Z);
[E,V] = infer(arima_model,Z);
h1 = kstest(E)
h2 = lbqtest(E)
parcorr(Z,'NumLags',300)

xxx=randn(10000,1);
kstest(xxx)

