function [Nh_SF,mse] = select_NhSF(EEG_raw, slicesPerSubts, a, trigs, EEG_input_data_ch, y_scsa)

% author: Gabriel Benigno

% select Nh_SF (ie best Nh value for shcrodinger filtering) by minimizing
% mean squared error Delta(Nh) (equations 7 of paper)

%%% inputs:
% EEG_raw: raw dataset
% slicesPerSubts: number of slice epochs per sub-timeseries
% a: index of slice epoch corresponding to start of the sub-timeseries. with slicesPerSubts, creates the sub-timeseries
% trigs: vector of slice acquisition trigger timings
% EEG_input_data_ch: input data for a single channel
% y_scsa: full reconstruction of signal

%%% outputs:
% Nh_SF: best Nh value for shcrodinger filtering
% mse: mean squared error Delta(Nh) (equation 7a of paper)


thr=0.025;
[t, ~, ~, t_shift, Ton_shift, Toff_shift] = int_thresh(EEG_raw, thr, slicesPerSubts, a, trigs);

Ton_vec = Tvec(Ton_shift,t_shift);
Toff_vec = Tvec(Toff_shift,t_shift);

y_temp=EEG_input_data_ch(t);
yoff=Toff_vec(:).*y_temp(:);
y = y_scsa;
yon = Ton_vec(:).*y;
[~,soff]=fft_function(1/2048,length(t),yoff);
L_Nh = size(y_scsa,2);
F = ceil(length(t)/2+0.5);
son = zeros(F,L_Nh);

for Nh = 1:L_Nh
    [f,son(:,Nh)]=fft_function(1/2048,length(t),yon(:,Nh));
end

fi=find(f<800000);

mse = mean((son(fi,:) - soff(fi)).^2)';

Nh_SF = find( mse == min(mse) );

end

