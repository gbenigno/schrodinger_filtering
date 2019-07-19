function [mra,mfsf] = run_perormance_evaluation(EEG_raw,steps,EEGins,trigs)

% author: Gabriel Benigno

% generate performance metrics ie meadian residual activity (MRA) and median fraction at slice frequencies (MFSF)

%%% inputs:
% EEG_raw: raw EEG dataset
% trigs: vector of slice acquisition trigger timings

%%% first generate steps and EEGins arrays eg:
% steps = ["aas","ds","sf","fastr","ica"];
% EEGins = {EEG_fmrib_aas EEG_fmrib_aas_spikeless ...
%     EEG_fmrib_sf_uninterp ...
%     EEG_fmrib_fastr ...
%     EEG_fmrib_aas_lpf150Hz_ica};
%
% another eg:
%
% steps = ["aas","sfwods","sfwds"];
% EEGins = {EEG_fmrib_aas ...
%     EEG_fmrib_sf_nods ...
%     EEG_fmrib_sf_uninterp};

%%% ouputs:
% mra: MRA value for each step
% rasf: MFSF value for each step

nch=30;

bands = [1 4; 4 8; 8 12; 12 30; 30 120];
band_names = ["theta", "delta", "alpha", "beta", "gamma"];

Tscan = length(trigs(1):trigs(end)); % number of time points during scanning

yon=zeros(nch,Tscan,length(steps));
yoff=zeros(size(yon));
son=zeros(nch,ceil(Tscan/2),length(steps));
soff=zeros(size(son));

s_band_mean = zeros(nch, length(band_names), length(steps) ); % ch, band, step
mra = NaN(size(s_band_mean)); % ch, band, step
f0=7;
n=16;
asf = zeros(nch, n+1, length(steps)); % ch, freq, step
mfsf = NaN(size(asf));

for step = 1:length(steps)
    
    [yon(:,:,step), son(:,:,step), yoff(:,:,step), soff(:,:,step), f_hires] = censored_spectrum(EEG_raw, EEGins{step}, trigs);
    s_band_mean(:,:,step) = sig_pres(bands,band_names,f_hires,soff(:,:,step));
    [fsl_vec, asf(:,:,step)] = art_rem(f0, n, f_hires, son(:,:,step));
    
    if step > 1
        mra(:,:,step) = median_residual_activity(s_band_mean(:,:,1), s_band_mean(:,:,step), band_names);
        mfsf(:,:,step) = slice_freq_ratio(fsl_vec, asf(:,:,1), asf(:,:,step));
    end

end

end

