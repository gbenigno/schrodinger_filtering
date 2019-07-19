function mra = median_residual_activity(bandstruc_ref,bandstruc,band_names)

% author: Gabriel Benigno

% calculate mra

%%% inputs
% bandstruc_ref: mean spectrum mag from sig_pres.m of reference eg aas
% bandstruc: eg sf
% names: see sig_pres.m

%%% outputs
% mra: rows=chs; columns=bands

mra = zeros(30, length(band_names) );

for band = 1:length(band_names)
    for ch = 1:30
        bandstruc_current = bandstruc(ch,band);
        bandstruc_ref_current = bandstruc_ref(ch,band);
        mra(ch,band) = (bandstruc_current - bandstruc_ref_current) / bandstruc_ref_current;
    end
end

end
