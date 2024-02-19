function [l,g,h] = LLKLIEP(theta,kP,kQ)
l = -mean(theta'*kP,2) + log(mean(exp(theta'*kQ),2));
% l1=-mean(theta'*kP,2);
% l2=log(sum(exp(theta'*kQ),2));%loss function
% l = l1+l2;
N_q = sum(exp(theta'*kQ),2);
g_q = exp(theta'*kQ)./ N_q;%RATIO
g = -mean(kP,2) + kQ*g_q';%gradient
    % hessian
if nargout>2
    HH = diag(g_q) - g_q'*g_q;
    h = kQ*HH*kQ';
end
end