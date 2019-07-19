function trigs = trig_vec(EEGin)

% author: Gabriel Benigno

% create vector of trigger timings

trigs_pre = EEGin.event;
trigs=zeros(length(trigs_pre) , 1);
for i = 1:length(trigs_pre)
    trigs(i)=trigs_pre(i).latency;
end

end

