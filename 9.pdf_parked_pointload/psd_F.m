function f = psd_F(X1,X2)
    % X1 = Dataset100(:,1,1);
    % X2 = Dataset200(:,1,1);
    nfft = 128;
    noverlap = 0;
    window = 256;%2^nextpow2(length(X1)/nfft);
    pxx1 = pwelch(X1,window,noverlap,nfft);
    pxx2 = pwelch(X2,window,noverlap,nfft);
%     m_pxx1 = mean(pxx1);
%     s_pxx1 = std(pxx1);
    % m_pxx2 = mean(pxx2);
    % s_pxx2 = std(pxx2);
    f = (sum(((pxx2-pxx1)).^2)/length(pxx1))/...
        (sum(((pxx1)).^2)/length(pxx2));
end