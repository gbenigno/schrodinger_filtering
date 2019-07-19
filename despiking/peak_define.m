function dataPks = peak_define(tq,realDataPkLocs,y_scsa_interp_dt,y_scsa_nn_interp)

% author: Gabriel Benigno

% define region of spike

%%% inputs:
% tq: time vector of y_scsa_interp_dt
% realDataPkLocs: indices of spikes
% y_scsa_interp_dt: interpolated and detrended version of y_scsa
% y_scsa_nn_interp: nnonegative and interpolated version of y_scsa

%%% outputs:
% dataPks: matrix of vectors; each vector is peak region padded by zeros





% step 1: find two intersections of base of spike with the trend
% step 2: create interval

dataPks = zeros(length(tq),length(realDataPkLocs));

for i = 1:length(realDataPkLocs)
    peak_idx_interp = find( abs(tq-realDataPkLocs(i)) == min(abs(tq-realDataPkLocs(i)))); % index of peak in interpolated data
    
    % find start of spike
    current = peak_idx_interp;
    while y_scsa_interp_dt(current) > 0
        current = current - 1;
    end
    start_of_peak = current;
    
    % find end of spike
    current = peak_idx_interp;
    while y_scsa_interp_dt(current) > 0
        current = current + 1;
    end
    end_of_peak = current;
    
    % create interval
    dataPks(start_of_peak:end_of_peak , i) = y_scsa_nn_interp(start_of_peak:end_of_peak);
    
end