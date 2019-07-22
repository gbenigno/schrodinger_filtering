# schrodinger_filtering

DATA AVAILABILITY

The raw dataset (referred to in the code as EEG_fmrib) is publicly available at https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/#install.

The post-AAS dataset (referred to in the code as EEG_fmrib_aas) was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, only the steps up to and including 'CalcAvgArt' were performed. See CleanEx1.m of the FACET repository for details.

The dataset (referred to in the code as EEG_fmrib_aas_pca_lpf150Hz_anc) with the full FASTR variant pipeline (i.e., AAS, OBS, low-pass filter at 150 Hz, and ANC) applied was obtained by processing EEG_fmrib using the FACET toolbox (https://github.com/hansiglaser/facet). In particular, all steps from the CleanEx1.m script in the FACET repository were performed, but applying the low-pass filter at 150 Hz rather than 70 Hz.



