function EEG_aas_spikeless = sf_ds(ch_vec,EEG_aas,EEG_raw,trigs)

% author: Gabriel Benigno

% perform de-spiking step of Schrodinger filtering

%%% inputs: 
% ch_vec: vector of channels to loop through
% EEG_aas: dataset following AAS
% EEG_raw: raw dataset
% trigs: vector of slice acquisition trigger timings

%%% output:
% EEG_aas_spikeless: output following de-spiking

EEG_aas_spikeless = EEG_aas;
EEG_aas_spikeless.data( : , (trigs(1) : trigs(end)) - 75 ) = NaN;

EEG_aas_spikeless_data = EEG_aas_spikeless.data;

slicesPerSubts=1; % for this step, only one slice epoch at a time is processed

for ch = ch_vec
    
    disp( ['ds ch ', num2str(ch)] )
    
    EEG_aas_spikeless_data_ch = EEG_aas_spikeless_data(ch,:);

    for subts = 1:840 % for each sub-ts (which in this case is a single epoch)

        [t, ~, ~, ~, Ton_shift, ~, y_aas] = subts_generate(EEG_raw, EEG_aas, ch, subts, slicesPerSubts, trigs);

        y_aas_spikeless_pos = sf_pos_despike(y_aas,Ton_shift); % remove positive spikes
        y_aas_spikeless_pos_flipped = -y_aas_spikeless_pos; % flip
        y_aas_spikeless_pos_flipped_spikeless_neg = sf_pos_despike(y_aas_spikeless_pos_flipped, Ton_shift); % remove negative spikes
        y_aas_spikeless_posneg = -y_aas_spikeless_pos_flipped_spikeless_neg; % flip back

        EEG_aas_spikeless_data_ch(t) = y_aas_spikeless_posneg;

    end
    
    EEG_aas_spikeless_data(ch,:) = EEG_aas_spikeless_data_ch;

end

EEG_aas_spikeless.data = EEG_aas_spikeless_data;

end

