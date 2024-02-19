function dist = acf_L2(X1,X2,k)

c1 = autocorr(X1,k);
c2 = autocorr(X2,k);
dist = norm(c1-c2);


end