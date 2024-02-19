function [KL_d] = density_ratio(xq,xp)
% xp = X1;xq = X2;
% xp = zscore(xp);xq = zscore(xq);
% xq = M_A2_slct; xp = M_A1_slct;
% xq = xq - mean(xp);xp = xp - mean(xp); %20191028 delete this line.
% xp = xp./std(xp); xq = xq./std(xq);
% xp = xp';xq = xq';
% dimension, np, nq
m = length(xp); [~,np] = size(xp); [~,nq] = size(xq);
kp = kernel_linear(xp); kq = kernel_linear(xq);

%Naive subgradient descent
theta = zeros(size(kq,1),1); lambda = 0*log(m)/sqrt(np);
% tic
step = 1; slength = inf; iter = 0; fold = inf;
while(slength > 1e-5)
    [f, gt] = LLKLIEP(theta,kp,kq); % return loss function and gradient
    g = zeros(size(gt)); 

    id = abs(theta)>0; % for theta>0 or theta<0
    g(id) = gt(id) + lambda*sign(theta(id)); 
    id = theta==0 & gt > lambda;% for theta==0
    g(id) = gt(id) - lambda;
    id = theta==0 & gt < -lambda;% for theta==0
    g(id) = gt(id) + lambda;
    theta = theta - step*g./(iter+1); % update theta
    slength = step*norm(g)./(iter+1);
    fdiff = abs(f - fold);

    %display some stuffs
    if iter > 50000
        disp('max iteration reached.')
        break;
    else
        iter = iter+1;
        fdiff = abs(f - fold);
        fold = f;
%         if ~mod(iter,100)
%             fprintf('%d, %.5f, %.5f, %.5f, nz: %d\n',...
%                 iter, slength,fdiff,full(fold),sum(theta(1:end-m)~=0))
%         end
    end
end
% toc
%% KL-divergence
[KL_d, ~,~] = LLKLIEP(theta,kp,kq); KL_d = -KL_d;% return loss function
% ra = mean(theta' * kp,2); 