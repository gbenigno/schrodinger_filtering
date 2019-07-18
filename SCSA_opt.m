function [y_scsa_opt, J_opt, Nh_vec_opt] = SCSA_opt(yin, h_vec, dx, gm)

% function to find mean squared error for different h-values

% inputs:
% yin: 1D input signal
% h_vec: 1D vector of h-values
% dx: time base (set to 1)
% gm: gamma (set to 0.5)

% outputs:
% y_scsa_opt: matrix of SCSA reconstructions of yin; one column for each h-value
% J_opt: vector of mean squared errors; one value for each h-value
% Nh_vec_opt: vector of Nh-values; one value for each h-value


J_opt = zeros(length(h_vec),1); % mean squared error
Nh_vec_opt = zeros(length(h_vec),1);
y_scsa_opt = zeros(length(yin) , length(h_vec));
M = length(yin); % make sure y is 1D before doing this
for i = 1:length(h_vec)
    
    [y_scsa_opt(:,i), Nh_vec_opt(i), ~, ~, ~] = SCSA(yin, h_vec(i), dx, gm);
    ytemp = y_scsa_opt(:,i);
    
    if size(ytemp) ~= size(yin)
        ytemp = ytemp';
    end
    J_opt(i) = (1/M)*sum( (yin - ytemp).^2);
end

end
