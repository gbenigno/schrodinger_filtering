function [EEG_fmrib,trigs,EEG_fmrib_aas,EEG_fmrib_fastr,EEG_fmrib_aas_lpf150Hz_ica] = sf_init(working_directory)

% author: Gabriel Benigno

% initializations: load datasets and generate vector of trigger timings

% input: absolute path to working directory that contains the datasets

% outputs: raw dataset (EEG_fmrib), trigger timing vector (trigs), dataset
% with AAS applied (EEG_fmrib_aas), dataset with full fastr pipeline
% applied (EEG_fmrib_fastr), dataset with full ICA pipeline applied
% (EEG_fmrib_aas_lpf150Hz_ica)



% load fmrib data (raw)
EEG_fmrib = load([working_directory,'EEG_fmrib.mat']);
EEG_fmrib = EEG_fmrib.('EEG_fmrib');
    
% high-pass-filter EEG_fmrib in freq domain at 1 Hz to remove baseline
% drifts
for ch=1:30
    o=EEG_fmrib.data(ch,:);
    o_filt = fftgausshp(o,1,2048);
    EEG_fmrib.data(ch,:) = o_filt;
end
    
% slice timing trigger vector
trigs = trig_vec(EEG_fmrib);
trigs = [trigs; trigs(end)+min(diff(trigs))];

% load fmrib data (processed up to aas)
EEG_fmrib_aas = load([working_directory,'/EEG_fmrib_aas.mat']);
EEG_fmrib_aas = EEG_fmrib_aas.('EEGout');

% load fmrib data (full fastr pipeline ie aas-->pca-->lpf(150Hz)-->anc)
EEG_fmrib_fastr = load([working_directory,'EEG_fmrib_aas_pca_lpf150Hz_anc.mat']);
EEG_fmrib_fastr = EEG_fmrib_fastr.('EEGout');

% load full ica pipeline (aas-->lpf(150Hz)-->ica)
EEG_fmrib_aas_lpf150Hz_ica=load([working_directory,'EEG_fmrib_aas_lpf150Hz_ica.mat']);
EEG_fmrib_aas_lpf150Hz_ica = EEG_fmrib_aas_lpf150Hz_ica.('EEG_fmrib_aas_lpf150Hz_ica');

end

