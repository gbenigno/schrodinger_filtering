% initializations
working_directory = '/Users/Gabriel/Dropbox/Grad_School_MBP/Schrodinger_Filtering/';
[EEG_fmrib,trigs,EEG_fmrib_aas,EEG_fmrib_fastr,EEG_fmrib_aas_lpf150Hz_ica] = sf_init(working_directory);

% %%%%%%%%%%%%%%%%%%% Schrodinger filtering %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    ch_vec=1:30;    
   
    % ds -> gf
    
    globfilt_w_ds = globfilt_struct_func(); % pre-allocation of globfilt
    
    seq_both = [1 1];
    [EEG_aas_spikeless, EEG_sf_with_ds, globfilt_w_ds] = sf_algo(ch_vec, EEG_fmrib_aas, [], EEG_fmrib, globfilt_w_ds, seq_both, trigs);
    
    % gf (no ds)
    
    globfilt_without_ds_Jul11 = globfilt_struct_func(); % pre-allocation of globfilt
    
    EEG_input_for_gf = EEG_fmrib_aas;
    
    seq_gfonly = [0 1];
    [~, EEG_sf_without_ds_Jul11, globfilt_without_ds_Jul11] = sf_algo(ch_vec,EEG_fmrib_aas,EEG_input_for_gf,EEG_fmrib,globfilt_without_ds_Jul11,seq_gfonly,trigs);


% %%%%%%%%%%%%%% PERFORMANCE EVALUATION %%%%%%%%%%%%%%%%%%%

steps = ["aas","ds","sf_w_ds","sf_wo_ds","fastr","ica"];
EEGins = {EEG_fmrib_aas ...
    EEG_fmrib_aas_spikeless ...
    EEG_sf_with_ds ...
    EEG_sf_without_ds ...
    EEG_fmrib_fastr ...
    EEG_fmrib_aas_lpf150Hz_ica};

[mra,mfsf] = run_perormance_evaluation(EEG_fmrib,steps,EEGins,trigs);
  