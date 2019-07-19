function mfsf = slice_freq_ratio(fsl_vec, asf_aas, asf_method)

% author: Gabriel Benigno

% ratio of spectral power density avg at slice freq+/- 1Hz for method over aas

%%% inputs:
% fsl_vec: vector of slice acqusition frequencies
% asf_aas: amplitude at slice frequencies of post-AAS signal; row=ch; col=freq
% asf_method: amplitude at slice frequencies of post-method signal; row=ch; col=freq

mfsf=zeros(30, length(fsl_vec));

for ch=1:30
    for i = 1:length(fsl_vec)
        mfsf(ch,i) = (asf_method(ch,i) / asf_aas(ch,i)).^2;
    end
end

end

