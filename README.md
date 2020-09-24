# Schrödinger filtering

**Schrödinger filtering** is a new signal processing technique based on **semi-classical signal analysis (SCSA)** (https://doi.org/10.1007/s00498-012-0091-1), which treats an input signal as an attractive potential in the one-dimensional semi-classical Schrödinger operator. The discrete spectrum of the operator is used for the signal's analysis. The basis functions of this eigenvalue decomposition are localized and pulse-shaped. Individually, these components naturally suit needs such as peak detection or removal. Together, they reconstruct the input with a preference for energy-dense input peaks. Based on the value of a single decomposition parameter, the user has control over how much of the energy-sparse parts of the input (*e.g.*, noise) are captured in the reconstruction. Schrödinger filtering is the term describing any filtering use of SCSA.

Here, Schrödinger filtering has been used for gradient artifact spike removal from **EEG** data acquired in the presence of **fMRI**. Schrödinger filtering handles the data after processing by **average artifact subtraction** (https://doi.org/10.1006/nimg.2000.0599), which contains residual artifact with time-domain spikes.

## Author

All code in this repository was written in **MATLAB** by **Gabriel Benigno** except where stated otherwise.

## Citing this work

Please cite the following article if you use this repository's code in any capacity:

...

## Data availability

The data used in the Schrödinger filtering paper comprises a freely available online dataset available at https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/#tutorial, a dataset for which permission is required (https://doi.org/10.1101/253047), and simulated data. The code used to generate the simulated data is within the `figures` folder.

## Contact

Feedback is welcome: gbenigno@uwo.ca.
