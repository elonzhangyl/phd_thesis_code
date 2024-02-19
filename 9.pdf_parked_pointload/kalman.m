clear;clc
t = (-2:0.001:3)'*5;
n = randn(size(t)); % 噪声
s = sin(t);         % 信号
x = s + n;          % 带噪信号
R = cov(n);         % 测量误差协方差

% 【1】预测误差比较大的时候
Q = 200;
y = KalmanFilter(x,Q,R);
e = s - y;
figure;
subplot(211);
plot(x,'color',[0.2 0.5 0.8],'linewidth',2);hold on;
plot(y,'color',[1 0.6 0],'linewidth',2);
plot(s,'color',[1 0.9 0],'linewidth',2);
legend('观测','滤波','真值','location','Best');
axis tight
subplot(212);
plot(e,'color',[0.2 0.5 0.8]);axis tight
legend('误差','location','Best');
axis tight

function T = KalmanFilter(s,Q,R)
% --------------------------------
%  T = KalmanFilter(s,Q,R)
%  s   ╬╗Êãú¼┴ð¤‗┴┐
%  Q   Î┤╠¼╬¾▓¯ð¡À¢▓¯¥Ï
%  R   ╣█▓Ô╬¾▓¯ð¡À¢▓¯¥Ï
% --------------------------------
%  ¢¿┴óÍ╩ÁÒÈ╦Â»À¢│╠:
%  S(k+1) = 1*S(k) + T*V(k) + 0.5*T^2*a
%  V(k+1) = 0*S(k) + 1*V(k) + T*a
%  ¢¿┴ó╣█▓ÔÀ¢│╠:
%  y(k+1) = S(k+1) + v(k+1)
%  ╝┤ú║
%  X(k+1) = A * X(k)   + G*w(k+1); Èñ▓Ô─úð═
%  y(k+1) = H * X(k+1) + v(k+1);   ╣█▓Ô─úð═


N = length(s);
T = 1;         %  ▓╔Ð¨╝õ©¶─¼╚¤╬¬1
A = [1 T;0 1]; %  Î┤╠¼Î¬Êã¥Ïı¾
G = [T^2/2;T]; %  ┐ÏÍã┴┐¥Ïı¾
H = [1 0];     %  ╣█▓Ô¥Ïı¾

% │§╩╝╗»Á┌Ê╗©÷Î┤╠¼
Xu = [s(1); 0];
Pu = [0 0;0 0];
I  = [1 0;0 1];
T  = zeros(N,1);

for i = 2:N
    Xp = A * Xu;
    Pp = A * Pu * A' + G * Q * G';
    K  = Pp * H' * ( H * Pp * H' + R)^-1;
    Xu = ( I - K * H ) * Xp + K * s(i);
    Pu = ( I - K * H ) * Pp;
    T(i) = Xu(1);
end
 
end
