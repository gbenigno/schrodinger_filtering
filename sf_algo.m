function [EEG_aas_spikeless,EEG_sf,globfilt] = sf_algo(ch_vec,EEG_aas,EEG_input_for_gf,EEG_raw,globfilt,seq,trigs)

% author: Gabriel Benigno

% performs Schrodinger filtering algorithm

%%% inputs:
% ch_vec: vector of channels to loop through
% EEG_aas: dataset following AAS
% EEG_input_for_gf: the dataset to be inputted for global filtering
% EEG_raw: raw dataset
% globfilt: structure array containing various data from global filtering step; input as blank
% seq: sequence of steps. ie [1 1] for both ds and gf; [1 0] for ds but not gf; and [0 1] for gf but not ds
% trigs: vector of slice acquisition trigger timings

%%% outputs:
% EEG_aas_spikeless: output following de-spiking
% EEG_sf: output following entire schrodinger filtering algorithm
% globfilt: structure array containing various data from global filtering step

    if seq(1)==1 % if despiking is selected

        %%%%%%%% STEP 1: DE-SPIKING %%%%%%%
        EEG_aas_spikeless = sf_ds(ch_vec,EEG_aas,EEG_raw,trigs);

        if seq(2)==1
            %% %%%%%% STEP 2: GLOBAL FILTERING (FOLLOWING DE-SPIKING) %%%%%%%
            EEG_input_for_gf = EEG_aas_spikeless;
            [EEG_sf,globfilt] = sf_gf(ch_vec,EEG_aas,EEG_input_for_gf,EEG_raw,globfilt,trigs);
        end

        % if seq(2)=0 within seq(1)=1, nothing is done

    elseif seq(1)==0 && seq(2)==1

        %%%%%%%%% ALTERNATE: GLOBAL FILTERING (WITHOUT DE-SPIKING) %%%%%%%
        [EEG_sf,globfilt] = sf_gf(ch_vec, EEG_aas, EEG_input_for_gf, EEG_raw, globfilt, trigs);
        EEG_aas_spikeless=[];

    end

end

