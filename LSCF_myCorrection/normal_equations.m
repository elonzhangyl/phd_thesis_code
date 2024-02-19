function [A, B] = normal_equations(P, T, n)
% [A, B] = normal_equations(P, T, n)
% solves normal equations of the least-squares problem
%
% P    matrix obtained from lscf.m
% T    matrix obtained from lscf.m
% n    at order n 
%
% A    coefficient of the numerator
% B    coefficient of the denominator

% normal equations
X = -real(P'*T) ;
Y = real(P'*P) ;
Z = real(T'*T) ;

% combinate normal equations
M = Z - X'*Y^(-1)*X ;

% linear system with B(end)=1 to avoid trivial solution
% nx1 vector
B = -M(1:n, 1:n) \ M(1:n, n+1) ;

% (n+1)x1 vector
B = [B; 1] ;

% (n+1)x1 vector
A = -Y\(X*B) ;

