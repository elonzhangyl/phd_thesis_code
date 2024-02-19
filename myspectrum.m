function [w, frf] = myspectrum(data, fs)
% [freq, fftdata] = spectrum(data, fs)
%
% data   time data vector
% fs     sample frequency in Hz

% ts = 1 / fs ;

L = length(data) ;

nfft = 2^nextpow2(L) ;

fftdata = fft(data, nfft) / L ; % normalization

fftdata = fftdata(1: nfft / 2) ;

% multiply by 2 to take into account the fact that we threw out second half
fftdata = fftdata * 2 ;
freq = fs / 2. * linspace(0., 1., nfft / 2) ;

w = 2*pi*freq ;
w = w.' ;
frf = fftdata ;

% returns natural frequency vector w in rad/s  and complex frequency response
% function frf

