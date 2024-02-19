%for decide AR order
clc;clear;

load('und_seed101To111.mat')
node = 1;
resam = 10;
dataset_ref = dataset1(1:resam:end,node,1);

Z = dataset_ref(:,node,1);
% Z = autocorr(Z,'NumLags',300-1);
% figure;
% subplot(2,1,1);autocorr(Z,'NumLags',size(Z,1)-1);
% subplot(2,1,2);parcorr(Z,'NumLags',100-1);
%estimate ARMA coefficients
p = 20;q = 1;

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

%% PLOT
xx=1:length(aic);yy=aic;
figure('position',[0 0 255 200],'units','points');
p = plot(xx,yy);
% p.Color = 'k';
p.LineStyle = '-';
p.LineWidth = 1;
p.Marker = 'o';
p.MarkerSize = 2;
% p.MarkerEdgeColor = 'k';
% p.MarkerFaceColor = 'k';
% xlim([1,250]);
% ylim([0 0.12]);
% xticks([1,50:50:300]);
xlabel('Model order');
ylabel('AIC');
% ytickformat('percentage')
% yline(0.05,'r');
set(findall(gcf,'-property','FontSize'),'FontSize',7)

