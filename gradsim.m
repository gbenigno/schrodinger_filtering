function [t, G, trigs, ppe, Ns] = gradsim(Fs, TRsl, Nsl, pausepercent, stf)
% simulate gradient artifact using sawtooth waves
% Gabriel Benigno; May 2020
% 
% INPUTS
% Fs - sampling frequency (Hz)
% TRsl - MRI slice repetition time
% Nsl - number of slices
% pausepercent - (e.g. 0.5 for 50%) % of each epoch that pauses
% stf - sawtooth frequency (Hz)
% 
% OUTPUTS
% tG - Time samples (s) of the simulated waveform. Begins at 0 s. 
% G - samples of gradient waveform (arbitrary unit) corresponding to tG
% trigs - vector of indices of the latencies of the triggers corresponding to the start of each slice acq

ppe = floor(TRsl*Fs); % points per epoch (incl pauses)
ppa = ceil((1-pausepercent)*ppe); ppa = ppa + mod(ppa,2) + mod(ppe,2); % points per k-space sampling per epoch
ppp = (ppe-ppa)/2; % points per pause; one at beginning and one at end of each epoch

t = 0 : 1/Fs : (ppe*Nsl-1)/Fs;
t = reshape(t',ppe,Nsl);

% logical vector corresponding to pausing between k-space acquisition
ip = false(size(t));
ip([1:ppp end-ppp+1:end] , :)=true; 

% logical vector corresponding to k-space acquisition samples
ia = ~ip;

% simulate gradient artifact separately for each epoch
G = zeros(size(t));

G(ia)=sawtooth(2*pi*stf*t(ia), 1/2);

% add "buffer" samples (samples before and after fMRI acquisition)
ppb = ppe; % points per buffer
bfr1=(0:ppb-1)'/Fs;
tmp=t;
t = [bfr1 tmp+ppb/Fs];
bfr2 = bfr1+numel(t)/Fs;
t = [t bfr2];

% polish tG and G arrays
t = t(:);
G=G(:); G = [zeros(ppb,1); G; zeros(ppb,1)];

% trigger time indices
trigs = (1 : ppe : Nsl*ppe)' + ppb ;

Ns = length(t);

end