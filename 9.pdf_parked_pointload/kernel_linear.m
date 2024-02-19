function [k] = kernel_linear(x)
[d n] = size(x); k = [];
tt = ones(d); tt = triu(tt,1); tt=tt(:); idx = tt~=0;
    for i = 1:n
        t = x(:,i) * x(:,i)';
        t = triu(t,1); t = t(:);
        k = [k, t(idx)];
    end
% k = [k;x.^2];
% k = [k;x.^2;x];
k = [k];

end
%off diag element first, and then diag elements
%this function used to generate features
