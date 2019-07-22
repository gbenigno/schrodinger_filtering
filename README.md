# schrodinger_filtering

DATA AVAILABILITY

The raw dataset (referred to in the code as EEG_fmrib) is publicly available at https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/#install.

The post-AAS dataset (referred to in the code as EEG_fmrib_aas) was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, only the steps up to and including 'CalcAvgArt' were performed. See CleanEx1.m of the FACET repository for details.

The dataset (referred to in the code as EEG_fmrib_aas_pca_lpf150Hz_anc) with the full FASTR variant pipeline (i.e., AAS, OBS, low-pass filter at 150 Hz, and ANC) applied was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, all steps from the CleanEx1.m script in the FACET repository were performed, but applying the low-pass filter at 150 Hz rather than 70 Hz.

The dataset with the full ICA pipeline applied to it is found at https://www.dropbox.com/s/yle0urcdygeikfl/ica.mat?dl=0 as a MATLAB structure array called ica.mat. The fields of this structure array are y, A, W, x, Z, and xc. y is the pre-ICA, post-AAS dataset. A is the mixing matrix. W is the unmixing matrix. x is the matrix of independent components. Z is the diagonal matrix containing a diagonal entry of 1 for ICs marked artifact and 0 for ICs marked signal. xc is y following artifact IC subtraction. ICA was performed using FastICA on MATLAB (https://research.ics.aalto.fi/ica/fastica/).
