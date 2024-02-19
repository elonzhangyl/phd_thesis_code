function dist = PEDiv(P,Q) %pearson divergence
% P = X1; Q = X2;

mu_p = mean(P);mu_q = mean(Q);
sigma_p = cov(P);sigma_q = cov(Q);

dist = mean(mvnpdf(P,mu_q,sigma_q)./mvnpdf(P,mu_p,sigma_p))-1;

end
