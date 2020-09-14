# Schrödinger filtering

**Schrödinger filtering** is a new signal processing technique based on **semi-classical signal analysis (SCSA)** (https://doi.org/10.1007/s00498-012-0091-1), which treats an input signal as an attractive potential in the one-dimensional semi-classical Schrödinger operator. The discrete spectrum of the operator is used for the signal's analysis. The basis functions of this eigenvalue decomposition are localized and pulse-shaped. Individually, these components naturally suit needs such as peak detection or removal. Together, they reconstruct the input with a preference for energy-dense input peaks. Based on the value of a single decomposition parameter, the user has control over how much of the energy-sparse parts of the input (*e.g.*, noise) are captured in the reconstruction. Schrödinger filtering is the term describing any filtering use of SCSA.

Here, Schrödinger filtering has been used for gradient artifact spike removal from **EEG** data acquired in the presence of **fMRI**. Schrödinger filtering handles the data after processing by **average artifact subtraction** (https://doi.org/10.1006/nimg.2000.0599), which contains residual artifact with time-domain spikes.

## Author

All code in this repository was written in **MATLAB** by **Gabriel Benigno** except where stated otherwise.

## Citing this work

Please cite the following article if you use this repository's code in any capacity:

...

## Data availability

The data used in the Schrödinger filtering paper comprised a freely available online dataset available at https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/#tutorial, a dataset for which permission is required (https://doi.org/10.1101/253047), and simulated data. The freely available real data are available in the `realdata` folder in csv files, which can be loaded in MATLAB using the `readmatrix` function. The first line of each csv file has details. This first line does not get loaded into MATLAB. The code used to generate the simulated data is provided in the `simulations` folder.

## Future improvements

1. This code was written with universality in mind. That is, it attempts to support a large subset of all possible datasets. If, as it is used on an increasing number of datasets, a need for greater generality becomes apparent, do not hesitate to contact the author (gbenigno@uwo.ca), and he will swiftly update the code accordingly. Also, feel free to collaborate by forking the repo.

2. There are plans to convert all toolbox functions into the C language, beginning with the least efficient ones, while keeping the MATLAB wrapper for convenience. Thus, near-future improvements on this front can be expected. Feel free to contribute to this process through GitHub, or to make specific requests to the author (gbenigno@uwo.ca).
