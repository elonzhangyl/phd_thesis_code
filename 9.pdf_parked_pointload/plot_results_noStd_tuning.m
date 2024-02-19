clear
load('result_pdf_KLDiv_tuning_damage2and4.mat')
a2 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1:5);
a4 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,6:10);
load('result_pdf_KLDiv_tuning_0.01.mat')
a1 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD;
load('result_pdf_KLDiv_tuning_0.03.mat')
a3 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD;
b11 = cat(3,a1(:,:,1:25,:,:,:,1),a1(:,:,1:25,:,:,:,2),a1(:,:,1:25,:,:,:,3),a1(:,:,1:25,:,:,:,4),a1(:,:,1:25,:,:,:,5));
b21 = cat(3,a2(:,:,1:25,:,:,:,1),a2(:,:,1:25,:,:,:,2),a2(:,:,1:25,:,:,:,3),a2(:,:,1:25,:,:,:,4),a2(:,:,1:25,:,:,:,5));
b31 = cat(3,a3(:,:,1:25,:,:,:,1),a3(:,:,1:25,:,:,:,2),a3(:,:,1:25,:,:,:,3),a3(:,:,1:25,:,:,:,4),a3(:,:,1:25,:,:,:,5));
b41 = cat(3,a4(:,:,1:25,:,:,:,1),a4(:,:,1:25,:,:,:,2),a4(:,:,1:25,:,:,:,3),a4(:,:,1:25,:,:,:,4),a4(:,:,1:25,:,:,:,5));
b12 = cat(3,a1(:,:,26:50,:,:,:,1),a1(:,:,26:50,:,:,:,2),a1(:,:,26:50,:,:,:,3),a1(:,:,26:50,:,:,:,4),a1(:,:,26:50,:,:,:,5));
b22 = cat(3,a2(:,:,26:50,:,:,:,1),a2(:,:,26:50,:,:,:,2),a2(:,:,26:50,:,:,:,3),a2(:,:,26:50,:,:,:,4),a2(:,:,26:50,:,:,:,5));
b32 = cat(3,a3(:,:,26:50,:,:,:,1),a3(:,:,26:50,:,:,:,2),a3(:,:,26:50,:,:,:,3),a3(:,:,26:50,:,:,:,4),a3(:,:,26:50,:,:,:,5));
b42 = cat(3,a4(:,:,26:50,:,:,:,1),a4(:,:,26:50,:,:,:,2),a4(:,:,26:50,:,:,:,3),a4(:,:,26:50,:,:,:,4),a4(:,:,26:50,:,:,:,5));
b1 = cat(3,b11(:,:,1:100,:,:,:),b12(:,:,1:100,:,:,:));
b2 = cat(3,b21(:,:,1:100,:,:,:),b22(:,:,1:100,:,:,:));
b3 = cat(3,b31(:,:,1:100,:,:,:),b32(:,:,1:100,:,:,:));
b4 = cat(3,b41(:,:,1:100,:,:,:),b42(:,:,1:100,:,:,:));
b = cat(7,b1,b2,b3,b4);

DD=4;NRT=200;SS=1;SN=1;NL=1;SQ=51;
SCORE_avg = mean(b,2);
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,SN,NL,SQ,DD);
for m = 1:SS
    for i = 1:SN
        for j = 1:NL
            for k = 1:SQ
                for dd = 1:DD
                    [~,~,~,auc] = perfcurve(labels,...
                        squeeze(SCORE_avg(m,1,:,i,j,k,dd)),1);
                    AUC(m,i,j,k,dd) = auc;
                end
            end
        end
    end
end


%% plot Dimensionality-ROC
figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.45]*36.36);
% axes('position',[0.1 0.1 0.9 0.9]);
xx = [1 5:5:250];
yy1 = squeeze(AUC(:,:,:,:,4));
yy2 = squeeze(AUC(:,:,:,:,3));
yy3 = squeeze(AUC(:,:,:,:,2));
yy4 = squeeze(AUC(:,:,:,:,1));
lw = 1;ms = 2.5;
plot(xx,yy1,'LineWidth',lw,...
    'Marker','s','MarkerSize',ms);hold on
% plot(xx,yy2,'LineWidth',lw,...
%     'Marker','s','MarkerSize',ms);hold on
plot(xx,yy3,'LineWidth',lw,...
    'Marker','v','MarkerSize',ms);hold on
% plot(xx,yy4,'LineWidth',lw,...
%     'Marker','v','MarkerSize',ms);hold off
xlim([1,200]);
ylim([0.48 1]);
xticks([1,50:50:300]);
xlabel('Dimensionality of PDF, {\it{\fontname{Times} k}}');
ylabel('ROC-AUC');
legend('4% stiffness reduction',...
    '2% stiffness reduction',...
    'location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'fig4.Dimensionality-ROC.eps','Resolution',1000)

