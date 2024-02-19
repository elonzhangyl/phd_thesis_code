clear
% extract result
SCORE = zeros(1,20,40,21,1,17,6);
for i = 1:21
    load(['results_KLEIP_PCA_tuning_sq_pc',...
        num2str(i+4),'.mat']);
    SCORE(:,:,:,i,:,:,:) = SCORE_WW_NRR_NRT_SN_NL_SQ_DD;
end
SCORE_WW_NRR_NRT_SN_NL_SQ_DD = SCORE(:,:,:,:,:,:,4);

wd_list = 10;%[0.7,0.8,0.9,1];
pca_order_list = 5:25;
sqenc_list = 40:10:200;
alpha_list = 0.032;
SQ = length(sqenc_list);
NL = length(alpha_list);
SN = length(pca_order_list);%seed
NRR = 20; % number of reference time series
NRT = 40; % number of test time series
WW = length(wd_list);
DD = 6;
labels = (1:NRT)>(NRT/2);
avg = squeeze(mean(SCORE_WW_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,:),2));
auc = zeros(SN,SQ);
for sq=1:SQ
    for sn=1:SN
        [~,~,~,auc_temp] = perfcurve(labels,squeeze(avg(:,sn,sq)),1);
        auc(sn,sq) = auc_temp; 
    end
end
% plot(sqenc_list,auc,'-o')
% mean(auc,2)
%% plot 3D bar
% Create the 3D bar chart
figurewidth = 9; %cm
figure1 = figure('Position',[10 10 figurewidth figurewidth]*29.1);
axes1 = axes('Parent',figure1);
contour(auc,10)
% xlabel('{\it{\fontname{Times} k}}','Position',...
%     [7.512749266607898,-1.911281214356563,-0.174529365211274])
% ylabel('{\it{\fontname{Times} d}}','Position',...
%     [-2.560610598753044,11.61227746034291,-0.135806307355715])
% zlabel('ROC-AUC')
% 
% % Change the x and y axis tick labels
% set(axes1,'PlotBoxAspectRatio',[1.06327947832136 1.61803398874989 1],...
%     'XTick',[1 4 7 10 13],'XTickLabel',{'40','55','70','85','100'},...
%     'YTick',[1 6 11 16 21],'YTickLabel',{'5','10','15','20','25'});
% % ylim([0 22])
% % xlim([2 14])
% view(axes1,[212.450886629956 31.8511783816966]);
hold on

%% plot line


%%
dam = 0.36;
load('und_seed101To400_1e6.mat')
dataset1 = dataset1(120*20+1:end,1,1);
wd_list = 10;%[0.7,0.8,0.9,1];
pca_order_list = 5:25;
sqenc_list = 40:10:200;
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

auc_line = zeros(size(auc,2),1);
for i = 1:size(auc,2)
    auc_line(i) = auc(pp(i)-4,i);
end

% axes('position',[0.1 0.1 0.9 0.9]);
xx = sqenc_list';
xx = floor(xx-40)/10+1;
yy = pp-4;
zz = auc_line;
p3 = plot3(xx,yy,zz);
p3.Color = 'r';
p3.LineStyle = '-';
p3.LineWidth = 2;
p3.Marker = 's';
p3.MarkerSize = 3;
p3.MarkerEdgeColor = 'r';
% p.MarkerFaceColor = 'w';
% xlim([40,110]);
% ylim([0.5 0.9]);
% xticks(40:20:120);
% xlabel('Sequence length, {\it{\fontname{Times} d}}');
% ylabel('ROC-AUC');

set(findall(gcf,'-property','FontSize'),'FontSize',7)


exportgraphics(figure1,'fig.sq_pc_owt.eps','Resolution',1000)