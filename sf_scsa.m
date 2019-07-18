function [h_and_Nh, y_scsa] = sf_scsa(EEG_raw,slicesPerSubts,EEG_input,a,ch,trigs)

% author: Gabriel Benigno

% perform SCSA and generate data of fully reconstructed signal as well as
% all h > h* and all Nh < Nh* 

%%% inputs:
% EEG_raw: raw dataset
% slicesPerSubts: number of slice epochs per sub-timeseries
% EEG_input: input dataset to be reconstructed by SCSA
% a: index of slice epoch corresponding to start of the sub-timeseries. with slicesPerSubts, creates the sub-timeseries
% ch: current channel
% trigs: vector of slice acquisition trigger timings

%%% outputs:
% h_and_Nh: two-column matrix: all h > h* and Nh < Nh*
% y_scsa: full reconstruction of input signal




% for first sub-ts, input empty matrix [] for h_SF argument
% for subsequent sub-ts, input the h_SF value gotten for first sub-ts

dx=1;
gm=0.5;

%% define t (interval of sub-ts) and Ton and Toff.
thr=0.025;
[t, ~, ~, ~, ~, ~] = int_thresh(EEG_raw, thr, slicesPerSubts, a, trigs);
yin = EEG_input.data(ch,t);

%% find Nh* for first sub-ts
h_vec=linspace(4.5,7.5,30)';

[y_scsa_opt, J, Nh_vec] = SCSA_opt(yin, h_vec, dx, gm);
idx = find(J == min(J)); 
h_star = h_vec(idx);
[~, Nh_star, ~, ~, ~] = SCSA(yin, h_star, dx, gm);

%% while finding Nh*, some Nh-values were < Nh*. record these and their h-values, and record their y_scsa data
h_and_Nh = [h_vec(idx:end) Nh_vec(idx:end)];
y_scsa = y_scsa_opt(:, idx:end);

%% (first sub-ts) acquire more data down until a low Nh-value

Nh_temp = 10^9;
h=h_and_Nh(end,1);
h0=h;
reps=0;

while Nh_temp >  round(Nh_star/30)
    [ys_temp, Nh_temp, ~, ~, ~] = SCSA(EEG_input.data(ch,t), h, dx, gm); % run scsa with h; get Nh_temp and ys_temp
    if (~ismember( Nh_temp, h_and_Nh(:,2) )) && (Nh_temp<Nh_star) % if Nh_temp is not a member of the computed Nh<Nh* and is < Nh*...
        reps=0;
        h_and_Nh = [h_and_Nh; [h Nh_temp]]; % ...add it and its h-value
        y_scsa = [y_scsa ys_temp'];
    else
        reps=reps+1; % number of times the same Nh results from a higher h-value
    end

%     % text to display
%     ch_str = [ 'ch = ', num2str(ch) ];
%     ch_str_w_space = white_space(18,ch_str);
%     h_vals = [ 'h = ' , num2str(h) ];
%     h_vals_w_space = white_space(18,h_vals);
%     Nh_temp_vals = [ 'Nh = ' , num2str(Nh_temp) , ' (', num2str(round(Nh_star/30)), ')' ];
%     Nh_temp_vals_w_space = white_space(18,Nh_temp_vals);
%     reps_vals = [ 'reps = ', num2str(reps) ];
%     disp([ch_str_w_space, h_vals_w_space , Nh_temp_vals_w_space , reps_vals ])

    if reps <= 3
        h = h+0.05*h0; 
    elseif (reps > 3)  && (reps <= 6)
        h = h+0.1*h0; 
    elseif (reps > 6) && (reps <= 9)
        h = h+0.15*h0;
    else
        h = h+0.2*h0;
    end
end

temp1=h_and_Nh(:,2);
temp2 = h_and_Nh;
h_and_Nh = NaN( max( temp1 ) , 2);
h_and_Nh(temp1, 1) = temp2(:,1);
h_and_Nh(temp1, 2) = temp2(:,2);

temp3 = y_scsa;
y_scsa = NaN(length(t),max(temp1));
y_scsa(:,temp1) = temp3;

end

