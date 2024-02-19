% clear
n_ob = 5;% each subgroup has n_ob observations
n_sam = 6;% subgroup
A3 = 1.427;
n_start = 60;% which tidal level to start

%% estimate expected value with Bayessian linear regression
% run extract freq
pred_m_owt2 = bayessian_linear...
    (freq_train_owt1(:),freq_train_owt2(:),...
    [freq_train_owt1(1:n_start+n_ob*n_sam)';...
    freq_test_owt1_dam(n_start+n_ob*n_sam+1:end)']);
pred_m_owt3 = bayessian_linear...
    (freq_train_owt2(:),freq_train_owt3(:),freq_train_owt2(:));
pred_m_owt1 = bayessian_linear...
    (freq_train_owt3(:),freq_train_owt1(:),freq_train_owt3(:));
    
%% X-bar control charts
figurewidth = 13; %cm
f = figure('Position',[10 10 figurewidth figurewidth*1]*29.1);

for i = 1:3
    subplot(3,1,i)
    if i == 1
        freq_test_und = freq_train_owt2(:);
        freq_test_dam = freq_train_owt2(:);
        pred_m = pred_m_owt2;
    elseif i == 2
        freq_test_und = freq_train_owt3(:);
        freq_test_dam = freq_train_owt3(:);
        pred_m = pred_m_owt3;
    elseif i == 3
        freq_test_und = freq_train_owt1(:);
        freq_test_dam = freq_test_owt1_dam(:);
        pred_m = pred_m_owt1;
    end
    [xbar,ucl,lcl,xbar_future] = xbar_analysis...
        (freq_test_und,freq_test_dam,pred_m,n_ob,n_sam,A3,n_start);
    %% plot
    xx = 1:2*n_sam;
    yy = [xbar,xbar_future];
    p2 = plot(xx(end/2:end),yy(end/2:end));hold on
    
    p2.LineStyle = '-';
    p2.LineWidth = 1;
    p2.Marker = 'o';
    p2.MarkerSize = 3;
    if i == 2
        p2.Color = '#00B050';
     p2.MarkerEdgeColor = '#00B050';
    else 
        p2.Color = '#FF0000';
        p2.MarkerEdgeColor = '#FF0000';
    end
    p2.MarkerFaceColor = 'w';
    p1 = plot(xx(1:end/2),yy(1:end/2));hold off
    p1.Color = '#00B050';
    p1.LineStyle = '-';
    p1.LineWidth = 1;
    p1.Marker = 'o';
    p1.MarkerSize = 3;
    p1.MarkerEdgeColor = '#00B050';
    p1.MarkerFaceColor = 'w';

    ylim([min(lcl-0.01,min(yy)-0.005),max(ucl+0.01,max(yy)+0.005)])
    yline(ucl,'r');
    yline(lcl,'r');
    xline(n_sam,'--')
    txt1 = sprintf('%.3f', ucl);
    txt2 = sprintf('%.3f', lcl);
    text(1,ucl+0.005,['UCL = ',txt1])
    text(1,lcl-0.005,['LCL = ',txt2])
    xlabel('Sample number')
    ylabel('$\bar{x}$','interpreter','latex')
    ttl = {'(a)','(b)','(c)'};
    title(ttl(i))
    set(findall(gcf,'-property','FontSize'),'FontSize',7)
end
