function dist = L2Div(P,Q)
% P = M_A1_slct; Q = M_A2_slct;
% P=X1;Q=X2;
mu_p = mean(P);mu_q = mean(Q);
sigma_p = cov(P);sigma_q = cov(Q);

y_p = mvnpdf(P,mu_p,sigma_p);
y_q = mvnpdf(P,mu_q,sigma_q);

dist = sum((y_p-y_q).^2);
% plot(1:17983,(y_p-y_q).^2)
end
