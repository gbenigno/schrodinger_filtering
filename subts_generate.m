function [t, Ton, Toff, t_shift, Ton_shift, Toff_shift, y_out] = subts_generate(EEGraw, EEGin, ch, a, slicesPerSubts,trigs)

% author: Gabriel Benigno

% generate a vector of the sub-timeseries of interest, along with various time vectors.

%%% inputs:
% EEG_raw: raw dataset
% EEGin: input dataset of interest
% ch: current channel
% a: index of slice epoch corresponding to start of the sub-timeseries. with slicesPerSubts, creates the sub-timeseries
% slicesPerSubts: number of slice epochs per sub-timeseries
% trigs: vector of slice acquisition trigger timings

%%% outputs:
% t: sub-timeseries
% Ton: parts of t during which gradients were on
% Toff: parts of t during which gradients were off
% t_shift, Ton_shift, Toff_shift: versions of the three outputs that are shifted according to t(1) = 1

thr=0.025;
[t, Ton, Toff, t_shift, Ton_shift, Toff_shift] = int_thresh(EEGraw, thr, slicesPerSubts, a,trigs);
y_out = EEGin.data(ch,t);

end

