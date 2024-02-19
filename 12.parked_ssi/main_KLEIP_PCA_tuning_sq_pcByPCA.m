clc
clear
dam = 0.36;
load('und_seed101To400_1e6.mat')
dataset1 = dataset1(120*20+1:end,1,1);
wd_list = 10;%[0.7,0.8,0.9,1];
pca_order_list = 5:25;
sqenc_list = 130;
alpha_list = 0.032;
NL = length(alpha_list);
SN = length(pca_order_list);%seed
NRR = 20; % number of reference time series
NRT = 40; % number of test time series
WW = length(wd_list);
DD = 6;

dataset_ref = dataset1;
SQ = length(sqenc_list);
pp = zeros(SQ,1);
for sq = 1:SQ
    sqenc = sqenc_list(sq);
    Dataset100 = dataset_ref;
    X1 = zeros(size(Dataset100,1)-sqenc,sqenc);
    for j=1:(size(Dataset100,1)-sqenc)
        X1(j,:) = Dataset100(j:j+sqenc-1);
    end
    [coeff,score,latent,~,x] = pca(X1);
%     var(:,)
    cum_perc = zeros(sqenc,1);
    for p = 1:sqenc
        cum_perc(p) = sum(x(1:p))/sum(x);
    end
%     plot(cum_perc)
    [RowNrs,ColNrs] = find(cum_perc>0.7 & cum_perc<0.9);
    dx = diff(x);
    pp(sq) = find(dx == min(dx(RowNrs)))+1;
%     plot(x,'-o');
%     xline(min(RowNrs));xline(max(RowNrs));
end
figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.6]*36.36);
p = plot(x);hold on
p.Color = [0, 0.4470, 0.7410];
p.LineStyle = '-';
p.LineWidth = 0.5;
p.Marker = 's';
p.MarkerSize = 4;
p.MarkerEdgeColor = [0, 0.4470, 0.7410];
set(findall(gcf,'-property','FontSize'),'FontSize',7)
% p.MarkerFaceColor = 'w';
plot([min(RowNrs);min(RowNrs)],[0;x(min(RowNrs))],'-r','linewidth',1);
hold on
plot([max(RowNrs);max(RowNrs)],[0;x(max(RowNrs))],'-r','linewidth',1);
hold on
scatter(pp,x(pp),'o','r')
hold on
xlim([1 30])
xlabel('PC number');
ylabel('Variance');
set(findall(gcf,'-property','FontSize'),'FontSize',7)


annotation('textarrow',[0.47 0.55],[0.27 0.43],...
    'headwidth',3,'headlength',8);hold on



rectangle('Position',[100 0.95 100 0.05])
hold on
axes('position',[0.38 0.5 0.5 0.4]);
xx = 1:max(sqenc_list);
yy = x;
% yy = yy(21:41);
p = plot(xx,yy);hold on
p.Color = [0, 0.4470, 0.7410];
p.LineStyle = '-';
p.LineWidth = 0.5;
p.Marker = 's';
p.MarkerSize = 4;
p.MarkerEdgeColor = [0, 0.4470, 0.7410];
% p.MarkerFaceColor = 'w';
xlim([9,17]);
ylim([0.9 2.6]);
scatter(pp,x(pp),'o','r')
xticks([10 13 16]);
yticks([1 2])
% xlabel('Dimensionality of PDF, $k$','Interpreter', 'latex');
% ylabel('ROC-AUC');

set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig.pc-ROC_owt.eps','Resolution',1000)