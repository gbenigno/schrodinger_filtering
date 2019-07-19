function s_band_mean = sig_pres(bands,names,f,s)

% author: Gabriel Benigno

% generate mean band fft mags to facilitate measuring signal preservation

%%% inputs
% bands: eg [1 4; 4 8; 8 12; 12 30; 30 130] gives 1-4Hz, 4-8Hz, etc
% names: eg ['theta' 'delta' 'alpha' 'beta' 'gamma']
% f: frequency vector of fft
% s: fft spectrum

%%% outputs
% s_band_mean: matrix. rows=chs; columns=bands

s_band_mean = zeros(30,length(names));

for i = 1:length(names) % for each band
    fb = find(  ( f >= bands(i,1) )  &  ( f < bands(i,2) )  ); % frequencies within band
    sb = s(:,fb); % spectral values within band
    for ch = 1:30
        s_band_mean(ch,i) = mean(sb(:,ch));
    end
end

end

