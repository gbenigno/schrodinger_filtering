# schrodinger_filtering

Schrödinger filtering is a novel processing technique for removing gradient artifact from EEG data that was collected during fMRI. It is based on semi-classical signal analysis (SCSA) (https://doi.org/10.1007/s00498-012-0091-1), which decomposes a signal according to its energy using the discrete spectrum of the semi-classical Schrödinger operator. Schrödinger filtering consists of two steps: (1) de-spiking, in which data spikes of the input time-domain dataset y_aas (raw dataset processed by average artifact subtraction) corresponding to residual gradient-related spikes are removed using SCSA; and (2) global filtering, in which SCSA is used to remove signal components from the timeseries of global residual artifact as well as to denoise the signal.

All code for Schrodinger filtering was written in MATLAB.

Please reference the URL to this repository if you use its code in any capacity.


======DATA AVAILABILITY======

The data used (in MATLAB .mat format) is available at https://www.dropbox.com/sh/hb7ivwjtkn5rnml/AAA6Eo7s5CdZK0AZdyAIuuS6a?dl=0.

The raw dataset EEG_fmrib is publicly available at https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/#install.

The post-AAS dataset EEG_fmrib_aas was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, only the steps up to and including 'CalcAvgArt' were performed. See CleanEx1.m of the FACET repository for details.

The dataset (EEG_fmrib_aas_pca_lpf150Hz_anc) with the full FASTR variant pipeline (i.e., AAS, OBS, low-pass filter at 150 Hz, and ANC) applied was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, all steps from the CleanEx1.m script in the FACET repository were performed, but with the low-pass filter applied at 150 Hz rather than 70 Hz.

The dataset (EEG_fmrib_aas_lpf150Hz_ica) with the full ICA pipeline applied to it is also found as a MATLAB structure array called ica.mat. The fields of this structure array are y, A, W, x, Z, and xc. y is the pre-ICA, post-AAS dataset (EEG channels 1-30; excluding EMG and ECG channels). A is the mixing matrix. W is the unmixing matrix. x is the matrix of independent components. Z is the diagonal matrix containing a diagonal entry of 1 for ICs marked artifact and 0 for ICs marked signal. xc is y following artifact IC subtraction. ICA was performed using FastICA on MATLAB (https://research.ics.aalto.fi/ica/fastica/).

EEG_aas_spikeless is EEG_fmrib_aas following the de-spiking step of Schrödinger filtering.

EEG_sf_with_ds is EEG_aas_spikeless following the global filtering step of Schrodinger filtering. globfilt_w_ds is a structure array containing various metadata from the global filtering step. In particular, it contains the fields y_scsa, h_and_Nh, h_and_Nh_SF, and mse. y_scsa is a 30x1 cell array, with each cell corresponding to one EEG channel. Each cell is size MxNh*, where M is the number of time points in the first sub-timeseries and Nh* is the number of Schrodinger components for full input signal reconstruction. Each column of each cell is a SCSA reconstruction using an Nh-value equal to the index of the column. h_and_Nh is also a 30x1 cell array. Each cell is a Nh*x2 matrix, where the first column is the h-value and the second column is the corresponding Nh_value. h_and_Nh_SF is a 30x2 matrix, with each row representing one of the channels and columns 1 and 2 representing h_SF and the corresponding value of Nh_SF, respectively. mse is a 30x1 cell array with each cell being a vector of length Nh*. The values of each entry of these vectors are the mean squared error Delta(Nh) and the index is the value of Nh.

EEG_sf_without_ds is EEG_fmrib_aas following the global filtering step of Schrödinger filtering. That is, it skips de-spiking. globfilt_without_ds is the corresponding structure array containing the metadata as described by globfilt_w_ds above but pertaining to EEG_sf_without_ds.


======RUNNING CODE======

Start at sf_master_script.m. Make sure to update the working directory absolute path (variable called wd) and that the datasets are found at wd. Also make sure that the paths for all code are added.
