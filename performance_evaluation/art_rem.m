function [fsl_vec,asf] = art_rem(f0, n, f, s)

% author: Gabriel Benigno

% generate mean fft mag per artifact slice frequency to facilitate
% measuring artifact removal

%%% inputs
% f0: fundamental artifact frequency 
% n: number of harmonics wish to look at
% f: frequency vector of s
% s: fft mag spectrum (f x nch)

%%% outputs:
% fsl_vec: vector of slice frequencies
% asf: amplitude at slice frequencies; row=ch; col=freq

fsl_vec = f0:f0:f0*(n+1);
fsl_vec = fsl_vec(:);
asf = zeros(30, n+1);

for ch=1:30
    for i = 1:n+1
        fi=f0*(1+(i-1));
        fsl = find( (f >= fi-1) & (f <= fi+1) );
        asf(ch,i) = mean(s(ch,fsl));
    end
end

end

