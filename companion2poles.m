function [poles] = companion2poles(B, n, dt)
% poles = companion2poles(B, n, dt) solves eigenvalue problem
% using companion matrix
%
% B       denominator coefficient
% n       at order n
% dt      sampling period
%
% poles   complex poles

% discrete poles
MC = [zeros(n-1, 1), eye(n-1); -B(1:n)'] ;
valp = eig(MC) ;

% poles
poles = log(valp)/dt ;

% stable poles with negative real part
poles = poles(real(poles) < 0) ;

% from conjugate to poles with positive imaginary part
poles = poles(imag(poles) > 0) ;
