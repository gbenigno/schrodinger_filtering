function y_spikeless = despike(y_scsa_nn_interp, y_comps_interp, spikeInfo, yin)

% author: Gabriel Benigno

% perform weighted subtraction of primary and secondary components 

%%% inputs:
% y_scsa_nn_interp: nnonegative and interpolated version of y_scsa
% y_comps_interp: interpolated version of y_comps
% spikeInfo: structure array containing information on spikes; output from higher code
% yin: input signal to scsa; used to get minimum value to add 

%%% output:
% y_spikeless: yin following weighted subtraction of primary and secondary components


if size(y_scsa_nn_interp)~=size(spikeInfo.dataPks(:,1))
    y_scsa_nn_interp=y_scsa_nn_interp';
end
    
y_spikeless=y_scsa_nn_interp;
for i = 1:length(spikeInfo.promComps) % for each spike
    % replace spike data (including tails) with straight line connecting ends
    vals = find(spikeInfo.dataPks(:,i) ~= 0); % indices of spikeInfo.dataPks(:,i) where the spike actually is
    x1=vals(1);
    x2=vals(end);
    y1=y_scsa_nn_interp(x1);
    y2=y_scsa_nn_interp(x2);
    m = (y2-y1)/(x2-x1);
    strline = m*vals;
    shift = y1 - strline(1);
    strline = strline+shift;
    
    %% spike removal
    if spikeInfo.accComps(i) == 0 % no accessory comp
        X = [y_scsa_nn_interp(vals) y_comps_interp(vals,spikeInfo.promComps(i)) ]; % column of ones, column of the data spike, column of the comp spike
        Y = strline; % the straight line
        B = X\Y;
        tilde = X*B;
        y_spikeless(vals) = tilde;
    else % there is an accessory comp
        X = [y_scsa_nn_interp(vals) y_comps_interp(vals,spikeInfo.promComps(i)) y_comps_interp(vals,spikeInfo.accComps(i))]; % column of ones, column of the data spike, column of the prom comp spike, column of the acc comp
        Y = strline; % the straight line
        B = X\Y;
        tilde = X*B;
        y_spikeless(vals) = tilde;
    end
end

y_spikeless = y_spikeless + min(yin);

% bring back to original resolution
x=(1:length(y_scsa_nn_interp))';
xq = linspace(x(1),x(end),length(yin))';
y_spikeless = interp1(x, y_spikeless, xq,'spline');

end

