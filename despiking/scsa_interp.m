function [tq, y_interp] = scsa_interp(y)

% author: Gabriel Benigno

% cubic-spline interpolate

%%% inputs: 
% y: input signal

%%% outputs:
% y_interp: interpolated version of y
% tq: time vector corresponding to y_interp, given that time vector of y is 1:length(y)

t=(1:length(y))';
tq = linspace(t(1),t(end),length(t)*100)';
y_interp = interp1(t, y, tq,'spline');

end

