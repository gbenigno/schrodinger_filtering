function [y_on, s_on, y_off, s_off, f_hires] = censored_spectrum(EEG_raw, EEG_in, trigs)

% author: Gabriel Benigno

% get spectrum of signal during Ton or Toff

%%% inputs
% EEG_raw: raw EEG dataset
% EEG_in: EEG input structure array
% trigs: vector of trigger slice acquisition timings

%%% outputs:
% y_on: signal with zeros during Toff
% y_off: signal with zeros during Ton
% s_on: FFT magnitude spectrum of y_on
% s_off: FFT magnitude spectrum of y_off
% f_hires: frequency vector of FFT magnitude spectrum

a=1;
slicesPerSubts=840;

thr=0.025;
[t, ~, ~, t_shift, Ton_shift, Toff_shift] = int_thresh(EEG_raw, thr, slicesPerSubts, a, trigs);
Ton_vec = Tvec(Ton_shift,t_shift);
Toff_vec = Tvec(Toff_shift,t_shift);

y_on = zeros(30,length(t));
y_off = zeros(30,length(t));

s_on = zeros( 30, ceil(length(t)/2) );
s_off = zeros( 30, ceil(length(t)/2) );

for ch = 1:30
    y = EEG_in.data(ch,t);
    y_on(ch,:) = y(:).*Ton_vec(:);
    y_off(ch,:) = y(:).*Toff_vec(:);
    
    [f_hires,s_on(ch,:)] = fft_function( 1/2048, length(t), y_on(ch,:) );
    [~,s_off(ch,:)] = fft_function( 1/2048, length(t), y_off(ch,:) );
    
end

end

