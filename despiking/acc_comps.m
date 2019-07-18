function accComps = acc_comps(promCompMaxIdxs,y_comps)

% author: Gabriel Benigno

% identify secondary components


% j is which prom comp (1st, 2nd, etc) the acc comp corresponds to

accComps = zeros(size(promCompMaxIdxs));

y_comps_norm = y_comps./max(y_comps);

minPkDist=0;
minPkProm=0.05;
minPkWidth=0;
minPkHeight=0;

for i = 1:size(y_comps_norm,2)
    [pks,locs]=findpeaks(y_comps_norm(:,i),'MinPeakProminence',minPkProm,'MinPeakDistance',minPkDist,'MinPeakWidth',minPkWidth,'MinPeakHeight',minPkHeight); % find peaks in comp that surpass the min prominence  
    
    
    if (length(locs)==2) && (abs(pks(1)-pks(2)) < 0.6) && (abs(locs(1)-locs(2)) < size(y_comps_norm,1)/6 ) % if there are only two prominent peaks for the comp

        j = find( (promCompMaxIdxs > locs(1)) & (promCompMaxIdxs < locs(2)) );
        if ~isempty(j)
            accComps(j) = i;
        end
    end
end

end

