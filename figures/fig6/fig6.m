clearvars;clc

TRsl = 0.1; pausepercent = 0.2; stf = 700; nch = 1; Nsl = 0; fs=2048;
ntrials=1000;

sar_sf = nan(ntrials,1);
sar_at = nan(ntrials,1);
sar_mf = nan(ntrials,1);
sar_wd = nan(ntrials,1);

sf_record = nan(408,ntrials);
x_record = nan(408,ntrials);
spks_record = nan(408,ntrials);
x_spks_record = nan(408,ntrials);

wd_record = nan(408,ntrials);
mf_record = nan(408,ntrials);
at_record = nan(408,ntrials);

[t, ~, ~, ~, Ns] = gradsim(fs, TRsl, Nsl, pausepercent, stf);

%%

cd ../..

parfor ii = 1:ntrials
    
    rng(ii)
    
    fprintf(' %s (ii=%u) \n',datetime, ii)
    
    wno = dsp.ColoredNoise('white',Ns,nch);
    x = wno();
    spks = zeros(size(x));
    bfr = round(0.05*length(t));
    jj=1;
    while jj < 5
        loc = randi([bfr, length(t)-bfr+1]);
        spktmp = 30*rms(x)*rand*gaussmf(t, [1/(2048), t(loc)]);
        if ~isempty(intersect(find(spks>0.01*30*rms(x)), find(spktmp>0.01*30*rms(x))))
            continue
        end
        spks = spks + spktmp;
        jj=jj+1;
    end
    jj=1;
    while jj < 5
        loc = randi([bfr, length(t)-bfr+1]);
        spktmp = 30*rms(x)*rand*gaussmf(t, [1/(2048), t(loc)]);
        if ~isempty(intersect(find(abs(spks)>0.01*30*rms(x)), find(abs(spktmp)>0.01*30*rms(x))))
            continue
        end
        spks = spks - spktmp;
        jj=jj+1;
    end
    
    x_spk = x+spks;

    sf = schrodingerFiltering(x_spk);
    at = filloutliers(x_spk,'clip');
    mf = medfilt1(x_spk);
    wd = x_spk - wdenoise(x_spk);
    
    sar_sf(ii) = sum((sf - spks).^2) ./ sum((sf - x).^2);
    sar_at(ii) = sum((at - spks).^2) ./ sum((at - x).^2);
    sar_mf(ii) = sum((mf - spks).^2) ./ sum((mf - x).^2);
    sar_wd(ii) = sum((wd - spks).^2) ./ sum((wd - x).^2);    
    
    sf_record(:,ii) = sf;
    x_record(:,ii) = x;
    spks_record(:,ii) = spks;
    x_spks_record(:,ii) = x_spk;
    wd_record(:,ii) = wd;
    at_record(:,ii) = at;
    mf_record(:,ii) = mf;
    
end

sf=sf_record;
x = x_record;
spks = spks_record;
x_spks = x_spks_record;
wd=wd_record;
mf=mf_record;
at=at_record;

clearvars -except t x x_spks sf wd mf at Ns fs sar_sf sar_at sar_mf sar_wd
cd figures/fig6
save('vars.mat')

%%
clearvars;clc
load('vars.mat')

ratio = sar_sf./sar_wd;
idx = find(ratio==max(ratio));

t = 1/fs : 1/fs : 408/fs;

figure
tl=tiledlayout(2,3);

nexttile(1)
plot(t,x_spks(:,idx),'k',t,x(:,idx),'r')
set(gca,'fontsize',24)

nexttile(2)
plot(t,x_spks(:,idx),'k',t,sf(:,idx),'r')
set(gca,'fontsize',24)

nexttile(3)
plot(t, x_spks(:,idx),'k',t, wd(:,idx), 'r')
set(gca,'fontsize',24)

nexttile(4)
plot(t,x_spks(:,idx),'k',t,mf(:,idx),'r')
set(gca,'fontsize',24)

nexttile(5)
plot(t,x_spks(:,idx),'k',t,at(:,idx),'r')
set(gca,'fontsize',24)

nexttile(6)
plot(sar_wd,sar_sf,'o',...
    sar_at,sar_sf,'o',...
    sar_mf, sar_sf, 'o', ...
    [0, max([sar_wd; sar_at; sar_mf])], [0, max([sar_wd; sar_at; sar_mf])], '--k', ...
    [0, max([sar_wd; sar_at; sar_mf])], 3*[0, max([sar_wd; sar_at; sar_mf])],'--k')
ylabel('signal-to-artifact ratio of SF');xlabel('signal-to-artifact ratio of other technique')
legend('comparison to WD','comparison to AT','comparison to MF','fontsize',18)
xlim([0 max(sar_wd)*1.01]);ylim([0 max(sar_sf)*1.01]);
set(gca,'fontsize',24)

tl.TileSpacing='none';




%%
function [t, G, trigs, ppe, Ns] = gradsim(Fs, TRsl, Nsl, pausepercent, stf)
% simulate gradient artifact using sawtooth waves
% Gabriel Benigno; May 2020
% 
% INPUTS
% Fs - sampling frequency (Hz)
% TRsl - MRI slice repetition time
% Nsl - number of slices
% pausepercent - (e.g. 0.5 for 50%) % of each epoch that pauses
% stf - sawtooth frequency (Hz)
% 
% OUTPUTS
% tG - Time samples (s) of the simulated waveform. Begins at 0 s. 
% G - samples of gradient waveform (arbitrary unit) corresponding to tG
% trigs - vector of indices of the latencies of the triggers corresponding to the start of each slice acq

ppe = floor(TRsl*Fs); % points per epoch (incl pauses)
ppa = ceil((1-pausepercent)*ppe); ppa = ppa + mod(ppa,2) + mod(ppe,2); % points per k-space sampling per epoch
ppp = (ppe-ppa)/2; % points per pause; one at beginning and one at end of each epoch

t = 0 : 1/Fs : (ppe*Nsl-1)/Fs;
t = reshape(t',ppe,Nsl);

% logical vector corresponding to pausing between k-space acquisition
ip = false(size(t));
ip([1:ppp end-ppp+1:end] , :)=true; 

% logical vector corresponding to k-space acquisition samples
ia = ~ip;

% simulate gradient artifact separately for each epoch
G = zeros(size(t));

G(ia)=sawtooth(2*pi*stf*t(ia), 1/2);

% add "buffer" samples (samples before and after fMRI acquisition)
ppb = ppe; % points per buffer
bfr1=(0:ppb-1)'/Fs;
tmp=t;
t = [bfr1 tmp+ppb/Fs];
bfr2 = bfr1+numel(t)/Fs;
t = [t bfr2];

% polish tG and G arrays
t = t(:);
G=G(:); G = [zeros(ppb,1); G; zeros(ppb,1)];

% trigger time indices
trigs = (1 : ppe : Nsl*ppe)' + ppb ;

Ns = length(t);

end