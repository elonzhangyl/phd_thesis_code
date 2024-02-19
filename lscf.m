function [fn, xin, frfnum, FST] = lscf(w, frf, n)
% [fn, xin, frfnum, FST] = lscf(w, frf, n)
% Linear Square Complex Frequency estimator using discrete-time z-model
%    
%             k=n
%            ----- 
%            \         k  
%            .     A z
%            /      
%            -----
%             k=0                            j k w dt
% H(z) = -------------------------,    z = e
%             k=n
%            ----- 
%            \         k    
%            .     B z
%            /      
%            -----
%             k=0
%
% w       natural frequency vector in rad/s
% frf     complex frequency response function vector
% n       if n is a scalar, direct identification at order n
%         if n is a vector, identification using stabilization chart
%
% fn      eigen frequency in Hz
% xin     modal damping factor
% frfnum  numerical identified complex frequency response function
% FST     cell array with frequency of stable poles in frequency and damping
%         empty is n is a scalar

w0 = w ;
L = length(w) ;
f = w/(2*pi) ;
yfrf = 20*log10(abs(frf)) ;
f0 = min(f) ;
fend = max(f) ;
w = 2*pi*(f-f0) ;
dt = 1/(2*(fend-f0)) ;

% cell array to save frequency of stable poles in frequency and damping
% at each iteration
FST = {} ;

% figure(1)
% hold on
% plot(w0/(2*pi), 20*log10(abs(frf)), 'LineWidth', 6, ...
% 'Color', [0.8, 0.8, 0.8], 'DisplayName', 'Measured FRF')

if length(n) == 1
%     disp('-----------------------------------------')
%     disp(['Identification at order ', num2str(n)])
%     disp('-----------------------------------------')
    P = zeros(L, n+1) ;
    for k=1:n+1
        % matrix for normals equations
        P(:, k) = zdomain(w, dt, k-1) ;
    end
    T = frf.*P ;
    % solve normal equations
    [A, B] = normal_equations(P, T, n) ;
    % numerical identified FRF
    frfnum = (P*A)./(P*B) ;
    % solve eigenvalue problem using companion matrix
    poles = companion2poles(B, n, dt) ;
    % modal parameters
    wn = abs(poles) ;
    wn = wn+2*pi*f0 ;
    xin = -real(poles)./wn ;
    % reorganize modal parameters
    [wn, iwn] = sort(wn) ;
    xin = xin(iwn)' ;
    fn = wn/(2*pi) ;
        
elseif length(n) > 1
%     disp('-----------------------------------------')
%     disp(['Identification between order ', num2str(n(1)), ...
%      ' and ', num2str(n(end))])
%     disp('-----------------------------------------')
    n_min = min(n);
    n_max = max(n);
    istab = (max(yfrf)-min(yfrf))/(n_max-n_min+1);
    Ptot = zeros(L, n_max+1);
    for k=1:n_max+1
        % matrix for normals equations
        Ptot(:, k) = zdomain(w, dt, k-1) ;
    end
    Ttot = frf.*Ptot ;
    ip = 0 ;
    iFST = 1 ;
    for p = n_min:n_max
        fn = [] ;
        xin = [] ;
        ff = [] ;
        xixi = [] ;
        mathp = [] ;
        P = Ptot(:, 1:p+1) ;
        T = Ttot(:, 1:p+1) ;
        % solve normal equations
        [A, B] = normal_equations(P, T, p) ;
        frfnum = (P*A)./(P*B) ;
        % solve eigenvalue problem using companion matrix
        poles = companion2poles(B, p, dt) ;
        % modal parameters
        wp = abs(poles) ;
        wp = wp+2*pi*f0 ;
        xip = -real(poles)./wp ;
        % reorganize modal parameters
        [wp, iwp] = sort(wp) ;
        xip = xip(iwp) ;
        fp = wp/(2*pi) ;
        % stabilization chart
        if ip == 0
            fmin1 = fp ;
            ximin1 = xip ;
        elseif ip > 0
            [fn, xin] = stabchart(fp, xip, fmin1, ximin1, yfrf, ip, ...
            f, fn, xin, ff, xixi, mathp, istab, p) ;
            xin = xin' ;
            fmin1 = fp ;
            ximin1 = xip ;
            if ~isempty(fn)
                FST{iFST} = fn ;
                iFST = iFST+1 ;
            end
        end
        ip = ip+1 ;
    end
end

% plot(w0/(2*pi), 20*log10(abs(frfnum)), ...
% 'k--', 'LineWidth', 2, 'DisplayName', 'Discrete FRF using z-model')
% xlabel('Frequency (Hz)')
% ylabel('FRF (dB)')
% box on
% legend('location', 'southwest')
    
    
    

