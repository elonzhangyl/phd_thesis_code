clear
%% AUC
load('result_pdf_KLDiv_k150_damage2and4.mat')
a2 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1:5);
a4 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,6:10);
b21 = cat(3,a2(:,:,1:25,:,:,:,1),a2(:,:,1:25,:,:,:,2),a2(:,:,1:25,:,:,:,3),a2(:,:,1:25,:,:,:,4),a2(:,:,1:25,:,:,:,5));
b41 = cat(3,a4(:,:,1:25,:,:,:,1),a4(:,:,1:25,:,:,:,2),a4(:,:,1:25,:,:,:,3),a4(:,:,1:25,:,:,:,4),a4(:,:,1:25,:,:,:,5));
b22 = cat(3,a2(:,:,26:50,:,:,:,1),a2(:,:,26:50,:,:,:,2),a2(:,:,26:50,:,:,:,3),a2(:,:,26:50,:,:,:,4),a2(:,:,26:50,:,:,:,5));
b42 = cat(3,a4(:,:,26:50,:,:,:,1),a4(:,:,26:50,:,:,:,2),a4(:,:,26:50,:,:,:,3),a4(:,:,26:50,:,:,:,4),a4(:,:,26:50,:,:,:,5));
b2 = cat(3,b21(:,:,1:100,:,:,:),b22(:,:,1:100,:,:,:));
b4 = cat(3,b41(:,:,1:100,:,:,:),b42(:,:,1:100,:,:,:));
b_pdf = cat(7,b2,b4);
b_pdf = squeeze(mean(b_pdf,2));
NRT = 200;SS = 4;NL = 4;DD = 2;
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NL,DD);
for ss = 1:SS
    for nl = 1:NL
        for dd = 1:DD
            [~,~,~,auc] = perfcurve(labels,...
                squeeze(b_pdf(ss,:,nl,dd)),1);
            AUC(ss,nl,dd) = auc;
        end
    end
end
AUC_pdf = squeeze(AUC);


load('result_ar_mah_arOrder25_ARslct13.mat')
a2 = SCORE_SS_NRT_NS_NL_SQ_DD(:,:,:,:,:,1:5);
a4 = SCORE_SS_NRT_NS_NL_SQ_DD(:,:,:,:,:,6:10);
b21 = cat(2,a2(:,1:25,:,:,:,1),a2(:,1:25,:,:,:,2),a2(:,1:25,:,:,:,3),a2(:,1:25,:,:,:,4),a2(:,1:25,:,:,:,5));
b41 = cat(2,a4(:,1:25,:,:,:,1),a4(:,1:25,:,:,:,2),a4(:,1:25,:,:,:,3),a4(:,1:25,:,:,:,4),a4(:,1:25,:,:,:,5));
b22 = cat(2,a2(:,26:50,:,:,:,1),a2(:,26:50,:,:,:,2),a2(:,26:50,:,:,:,3),a2(:,26:50,:,:,:,4),a2(:,26:50,:,:,:,5));
b42 = cat(2,a4(:,26:50,:,:,:,1),a4(:,26:50,:,:,:,2),a4(:,26:50,:,:,:,3),a4(:,26:50,:,:,:,4),a4(:,26:50,:,:,:,5));
b2 = cat(2,b21(:,1:100,:,:),b22(:,1:100,:,:));
b4 = cat(2,b41(:,1:100,:,:),b42(:,1:100,:,:));
% b = cat(7,b1,b2,b3,b4);
b_ar = cat(5,b2,b4);
NRT = 200;SS = 4;NL = 4;DD = 2;
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,1,NL,DD);
for ss = 1:SS
    for nl = 1:NL
        for dd = 1:DD
            [~,~,~,auc] = perfcurve(labels,...
                squeeze(b_ar(ss,:,1,nl,dd)),1);
            AUC(ss,1,nl,dd) = auc;
        end
    end
end
AUC_ar = squeeze(AUC);

load('result_acfL2_order90.mat')
a2 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,1:5);
a4 = SCORE_SS_NRR_NRT_SN_NL_SQ_DD(:,:,:,:,:,:,6:10);
b21 = cat(3,a2(:,:,1:25,:,:,:,1),a2(:,:,1:25,:,:,:,2),a2(:,:,1:25,:,:,:,3),a2(:,:,1:25,:,:,:,4),a2(:,:,1:25,:,:,:,5));
b41 = cat(3,a4(:,:,1:25,:,:,:,1),a4(:,:,1:25,:,:,:,2),a4(:,:,1:25,:,:,:,3),a4(:,:,1:25,:,:,:,4),a4(:,:,1:25,:,:,:,5));
b22 = cat(3,a2(:,:,26:50,:,:,:,1),a2(:,:,26:50,:,:,:,2),a2(:,:,26:50,:,:,:,3),a2(:,:,26:50,:,:,:,4),a2(:,:,26:50,:,:,:,5));
b42 = cat(3,a4(:,:,26:50,:,:,:,1),a4(:,:,26:50,:,:,:,2),a4(:,:,26:50,:,:,:,3),a4(:,:,26:50,:,:,:,4),a4(:,:,26:50,:,:,:,5));
b2 = cat(3,b21(:,:,1:100,:,:,:),b22(:,:,1:100,:,:,:));
b4 = cat(3,b41(:,:,1:100,:,:,:),b42(:,:,1:100,:,:,:));
b_acf = cat(7,b2,b4);
b_acf = squeeze(mean(b_acf,2));
NRT = 200;SS = 4;NL = 4;DD = 2;
labels = (1:NRT)>(NRT/2);
AUC = zeros(SS,NL,DD);
for ss = 1:SS
    for nl = 1:NL
        for dd = 1:DD
            [~,~,~,auc] = perfcurve(labels,...
                squeeze(b_acf(ss,:,nl,dd)),1);
            AUC(ss,nl,dd) = auc;
        end
    end
end
AUC_acf = squeeze(AUC);
%% dash line plot
for dd=1:2
figurewidth = 9; %cm
f = figure('outerPosition',[10 10 figurewidth figurewidth*1.1]*29.1);
% subplot(2,2,1);subplot(2,2,2);
plot(max(AUC_pdf(:,:,dd)),...
    '-o','MarkerSize',5,'linewidth',1);hold on
plot(max(AUC_ar(:,:,dd)),...
    '-+','MarkerSize',5,'linewidth',1);hold on
plot(max(AUC_acf(:,:,dd)),...
    '-v','MarkerSize',5,'linewidth',1);hold off
    ylim([0.6,1])
%     yticks(0.5:0.1:1)
%     xlim([0.9,4.1])
    xticks(1:4)
    xticklabels({'0%','15%','30%','45%'})
    xlabel('Noise level')
    ylabel('ROC-AUC')
tt = ['a','b','c','d'];
legend1 = legend('PDF-based','AR-based','ACF-based',...
        'location','southwest','Orientation','vertical');
set(findall(gcf,'-property','FontSize'),'FontSize',7)
exportgraphics(f,['fig7',tt(dd),'.comparison_errorbar',...
    num2str(dd),'.eps'],'Resolution',1000)
end

%% table
n1 = squeeze(cat(2,AUC_pdf(:,:,1),AUC_pdf(:,:,2)));
n2 = squeeze(cat(2,AUC_ar(:,:,1),AUC_ar(:,:,2)));
n3 = squeeze(cat(2,AUC_acf(:,:,1),AUC_acf(:,:,2)));
n = cat(1,n1,n2,n3);


