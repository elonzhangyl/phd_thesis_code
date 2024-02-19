function dist = PDiv(P,Q)
% P = M_A1_slct; Q = M_A2_slct;
beta = 0.5;
% P=X1;Q=X2;
mu_p = mean(P);mu_q = mean(Q);
sigma_p = cov(P);sigma_q = cov(Q);

y_p = mvnpdf(P,mu_p,sigma_p);
y_q = mvnpdf(P,mu_q,sigma_q);

w = y_p./(beta*y_q+(1-beta)*y_p);

dist = 1/length(P)*sum(w)-1;
% plot(1:17983,w)
end
