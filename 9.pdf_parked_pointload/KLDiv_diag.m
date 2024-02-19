function dist = KLDiv_diag(P,Q)
% P = M_A1_slct; Q = M_A2_slct;
mu_p = mean(P);mu_q = mean(Q);
sigma_p = cov(P);sigma_q = cov(Q);

% diag_p = diag(X10);diag_q = diag(X20);

for i = 1:size(P,2)
    sigma_p(i,i) = sigma_q(i,i);
end
    

dist = 1*(...
    trace(inv(sigma_p)*sigma_q)...
    +(mu_p-mu_q)*inv(sigma_p)*(mu_p-mu_q)'...
    -size(P,2)...
    +log(det(sigma_p)/det(sigma_q))...
    );

end
