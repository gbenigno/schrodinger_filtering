function [promComps,promCompMaxIdxs, realDataPkLocs] = prom_comps(y_comps,data_pk_locs, Ton_shift)

% author: Gabriel Benigno

% identify the indices of components that have only one prominent peak within Ton_shift ie identify primary components

%%% Inputs: 
% y_comps: matrix of Schrodinger components of signal of interest
% data_pks: values of the potential spikes
% data_pk_locs: indices of the potential spikes
% Ton_shift: portions of y_aas during which gradients were on; re-indexed such that first data point of y_aas is index 1

%%% outputs:
% promComps: indices of primary components eg first component out of Nh components
% promCompMaxIdxs: indices of maximum of primary components eg 34th data point in the component vector
% realDataPkLocs: indices of spikes in the EEG signal

promComps = [];
promCompMaxIdxs = [];
realDataPkLocs=[];

box = zeros(size(y_comps(:,1),1),1);
box(Ton_shift) = ones(1,length(Ton_shift));

if size(box) ~= size(y_comps(:,1))
    box=box';
end

for i = 1:round(size(y_comps,2)/3)

    minPkDist=0;
    minPkProm=0.5;
    minPkWidth=0;
    minPkHeight=10;%mean(y);
    thr=0;
    
    % find peak locations of comps
    [~,locs]=findpeaks(y_comps(:,i).*box,'MinPeakProminence',minPkProm,'MinPeakDistance',minPkDist,'MinPeakWidth',minPkWidth,'MinPeakHeight',minPkHeight,'Threshold',thr);    

    if length(locs)==1 % if there is only one peak in the comp 
                
        differ=abs(locs - data_pk_locs); % difference between comp peak location and data peak locations
        real_pk = find(differ<=1); % find which data peak (1st, 2nd, etc) overlaps with current comp peak
        if length(real_pk)==1 % if this peak overlaps with the data peak
            
            minPkDist=0;
            minPkProm=0.5;
            minPkWidth=0;
            minPkHeight=0;
            thr=0;
            
            [~,locs2]=findpeaks(y_comps(:,i).*box,'MinPeakProminence',minPkProm,'MinPeakDistance',minPkDist,'MinPeakWidth',minPkWidth,'MinPeakHeight',minPkHeight,'Threshold',thr);
            
            if length(locs2)==1
                realDataPkLocs = [realDataPkLocs; data_pk_locs(real_pk)];
                promComps = [promComps; i]; % record which comp this is
                promCompMaxIdxs = [promCompMaxIdxs; locs];
            end
            

        end
    end

end

end

