function [KL_d] = density_ratio_fmin(xq,xp)
% xp = X1;xq = X2;
% xp = zscore(xp);xq = zscore(xq);
% xq = M_A2_slct; xp = M_A1_slct;
% xq = xq - mean(xp);xp = xp - mean(xp); %20191028 delete this line.
% xp = xp./std(xp); xq = xq./std(xq);
% xp = xp';xq = xq';
% dimension, np, nq
kp = kernel_linear(xp); kq = kernel_linear(xq);

%% some settings for optimizer
opt = optimset('fminunc');
% opt.GradObj = 'on';
opt.MaxFunEvals = 10000;
opt.MaxIter = 10000;
opt.Display = 'off';
theta0 = zeros(size(kq,1),1); 
loss = @(theta) -mean(theta'*kp,2) + log(mean(exp(theta'*kq),2));
[theta_hat,~] = fminunc(loss,theta0,opt);
% toc
%% KL-divergence
KL_d = -loss(theta_hat);% return loss function
% ra = mean(theta' * kp,2); 