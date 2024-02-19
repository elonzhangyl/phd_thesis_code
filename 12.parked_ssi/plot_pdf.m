clc
clear
dam = 0.48;
load('und_seed101To400_1e6.mat')
dataset1 = dataset1(120*20+1:120*20+20*9*60,1,1);

load(['dam_scour',num2str(dam),'m_seed201To300_1e6.mat'])
dataset2 = dataset2(120*20+1:120*20+20*9*60,1,1);

pca_order = 14;
sqenc = 130;%5:5:80;

% corrcoef(dataset_und1(1:end-10,1),dataset_und1(10+1:end,1));
data1 = dataset1(1:end,1);
data2 = dataset2(1:end,1);
% acf1 = autocorr(data1,'NumLags',100);
% acf2 = autocorr(data2,'NumLags',100);
% acf11 = acf1*var(data1);
% acf22 = acf2*var(data2);
% plot(1:length(acf11),acf11,1:length(acf22),acf22)

X1 = zeros(size(data1,1)-sqenc,sqenc);
for j=1:(size(data1,1)-sqenc)
    X1(j,:) = data1(j:j+sqenc-1);
end
[coeff,score,latent,~,explain] = pca(X1);
X1_PCA = score(:,1:pca_order);

X2 = zeros(size(data2,1)-sqenc,sqenc);
for j=1:(size(data2,1)-sqenc)
    X2(j,:) = data2(j:j+sqenc-1);
end
[coeff,score,latent,~,explain] = pca(X2);
X2_PCA = score(:,1:pca_order);

figurewidth = 4.5; %cm
f1 = figure('Position',[10 10 figurewidth figurewidth]*36.36);
x = zscore(X1_PCA(:,[1,2]));
% ksdensity(x,'PlotFcn','contour')
[f,xi] = ksdensity(x);
X = reshape(xi(:,1),30,30);
Y = reshape(xi(:,2),30,30);
Z = reshape(f,30,30);
s = surf(X,Y,Z);
xlim([-3,3])
ylim([-3,3])
zlim([0,0.15])
% colormap('gray')
xlabel('1st PC','Rotation',20,'position',[-1.5,-5.5,0])
ylabel('2nd PC','Rotation',-35,'position',[-5.5,-2.8,0])
 text('position',[-2.8,2.8,0.14],'Interpreter','latex',...
     'String','\it p{\boldmath$(x)$}')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f1,'fig.pdf_owt_und.eps','Resolution',1000)


figurewidth = 4.5; %cm
f2 = figure('Position',[10 10 figurewidth figurewidth]*36.36);
x = zscore(X2_PCA(:,[1,2]));
% ksdensity(x,'PlotFcn','contour')
[f,xi] = ksdensity(x);
X = reshape(xi(:,1),30,30);
Y = reshape(xi(:,2),30,30);
Z = reshape(f,30,30);
s = surf(X,Y,Z);
xlim([-3,3])
ylim([-3,3])
zlim([0,0.15])
% colormap('gray')
xlabel('1st PC','Rotation',20,'position',[-1.5,-5.5,0])
ylabel('2nd PC','Rotation',-35,'position',[-5.5,-2.8,0])
 text('position',[-2.8,2.8,0.14],'Interpreter','latex',...
     'String','\it p{\boldmath$(x)$}')
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f2,'fig.pdf_owt_dam.eps','Resolution',1000)
