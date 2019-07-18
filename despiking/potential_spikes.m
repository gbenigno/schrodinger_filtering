function [data_pks,data_pk_locs] = potential_spikes(y_scsa,Ton_shift)

% author: Gabriel Benigno

% identify potential spikes in y_scsa

%%% inputs:
% y_scsa: input signal
% Ton_shift: % Ton_shift: portions of y_aas during which gradients were on; re-indexed such that first data point of y_aas is index 1

%%% outputs:
% data_pks: values of the potential spikes
% data_pk_locs: indices of the potential spikes

%%

warning('off','signal:findpeaks:largeMinPeakHeight') % turn off warnings of invalid min peak height

MAD = mad(y_scsa);

minPkDist=0;
minPkProm=20;
minPkWidth=0;
minPkHeight=3.5*MAD;
thr=0;

box = zeros(length(y_scsa),1);
box(Ton_shift(1:30)) = 1;
box(Ton_shift(end-29:end)) = 1;

if size(box) ~= size(y_scsa)
    box=box';
end

temp = y_scsa.*box;

[data_pks,data_pk_locs]=findpeaks(temp,'MinPeakProminence',minPkProm,'MinPeakDistance',minPkDist,'MinPeakWidth',minPkWidth,'MinPeakHeight',minPkHeight,'Threshold',thr);

data_pks=data_pks';
data_pk_locs=data_pk_locs';

end

