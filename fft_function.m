function [f,P1] = fft_function(TR,L,X)

% function to create real part of FFT magnitude spectrum and accompanying
% frequency vector in time units

%%% inputs:
% TR: sampling interval
% L: number of data points
% X: vector of data

%%% outputs:
% f: frequency vector with units of inverse of units of TR
% P1: real part of FFT magnitude spectrum of X

if mod(L,2)~=0
    X(end) = [];  
    L=L-1;
end

Fs=1/TR;
t = (0:L-1)*TR;

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
end

