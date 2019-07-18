function [t, Ton, Toff, t_shift, Ton_shift, Toff_shift] = int_thresh(EEGin, thr, slicesPerSubts, a, trigs)

% author: Gabriel Benigno

% using intensity thresholding, determine Ton and Toff

%%% inputs: 
%     EEGin: 
% EEGlab-format array that includes raw data (before any
% processing) and slice triggers. 
% EEGin.data should be size mxn, where m is # of channels and n is # of time points.
% Slice triggers are assumed to be correct. and under sub-label 'latency'
%     thr:
% threshold (% of max amplitude e.g. 0.035)
%     slicesPerSubts:
% number of slice epochs per sub-timeseries. 6 is a good number
%     a:
% index of starting slice epoch of sub-timeseries

%%% outputs
%    t:
% sub-timeseries
%   Ton and Toff: 
% non-contiguous intervals during which the gradients are on and off,
% respectively.

%% define input and sub-timeseries using slice triggers
y=EEGin.data(1,:);
b=a+slicesPerSubts;
t=trigs(a) : trigs(b);
t = t - 75;
t_shift = t - t(1) + 1;
t=t(:);

%% de-mean and normalize relative to max value
y_demean=detrend(y(t));
y_demean=y_demean./max(y_demean);
y_demean=y_demean(:);

%% initial Ton value by thresholding
Ton_shift = find( abs(y_demean) > thr );
Ton = Ton_shift + t(1) - 1;

%%
ends=12;
if slicesPerSubts == 1
  Ton = (Ton(1) : Ton(end));
  Ton1 = Ton(1);
  Tonend=Ton(end);
  Ton = [ (Ton1-ends):(Ton1-1) Ton (Tonend+1):(Tonend+ends) ]; % add values to ends
  Ton=Ton(:);
  Ton_shift = Ton - t(1) + 1;
else
    dTon=diff(Ton); % first diffs
    ends_of_artifact = find( dTon > 20 ); % indices of Ton for which the following value was more than 20 greater; means end of artifact (does not include final artifact)
    arts = cell(slicesPerSubts,1);
    
    arts{1} = Ton(1):Ton(ends_of_artifact(1));
    arts{1} = [(arts{1}(1)-ends):(arts{1}(1)-1) arts{1} (arts{1}(end)+1):(arts{1}(end)+ends) ];
    
    for i = 2:slicesPerSubts-1
        arts{i} = Ton(ends_of_artifact(i-1)+1) : Ton(ends_of_artifact(i));
        arts{i} = [(arts{i}(1)-ends):(arts{i}(1)-1) arts{i} (arts{i}(end)+1):(arts{i}(end)+ends) ];
    end
    arts{slicesPerSubts} = Ton(ends_of_artifact(slicesPerSubts-1)+1) : Ton(end);
    arts{slicesPerSubts} = [(arts{slicesPerSubts}(1)-ends):(arts{slicesPerSubts}(1)-1) arts{slicesPerSubts} (arts{slicesPerSubts}(end)+1):(arts{slicesPerSubts}(end)+ends) ];
    
    arts=arts';
    Ton=sort(cell2mat(arts));
    Ton_shift = Ton - t(1) + 1;
    
end

Toff = setdiff(t,Ton);
Toff_shift = Toff - t(1) + 1;

end

