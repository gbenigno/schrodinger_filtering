function [trend,y_dt] = scsa_detrend(y,opol)

% author: Gabriel Benigno

% de-trend input signal by a polynomial of specified order

%%% inputs: 
% y: input signal
% opol: polynomial order of trend

%%% outputs:
% trend: vector of trend
% y_dt: de-trended version of y

t=1:length(y);
if size(t) ~= size(y)
    t=t';
end
[p,~,mu] = polyfit(t,y,opol);
trend = polyval(p,t,[],mu);
y_dt = y - trend;

end

