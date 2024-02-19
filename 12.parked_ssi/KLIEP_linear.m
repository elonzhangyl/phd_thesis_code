function wh_x = KLIEP_linear(xq,xp)
% xp = x_nu;xq = x_de;
m = size(xp,1); np = size(xp,2); nq = np;
kp = kernel_linear(xp); kq = kernel_linear(xq);
%Naive subgradient descent
theta = sparse(zeros(size(kq,1),1)); lambda = 0.5*log(m)/sqrt(np);
step = 1; slength = inf; iter = 0; fold = inf;
ff=[];
% tic
while(slength > 1e-5)
%     iter
    [f, gt] = LLKLIEP(theta,kp,kq); % return loss function and gradient
    ff = [ff,f];
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
    if iter > 1000
        disp('max iteration reached.')
        break;
    else
        iter = iter+1;
        fdiff = abs(f - fold);
        fold = f;
%         if ~mod(iter,10)
%             fprintf('%d, %.5f, %.5f, %.5f, nz: %d\n',...
%                 iter, slength,fdiff,full(fold),sum(theta(1:end-m)~=0))
%         end
    end
end
% myt = toc
wh_x = -f;