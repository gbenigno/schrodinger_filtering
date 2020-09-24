clearvars;clc
load('../../realdata/data1/raw/ch001.mat')
load('../../realdata/data1/raw/ch007.mat')

%%
tt = 58887 : 304369+300;
tmp = double(ch007(tt));
TF = isoutlier(tmp);
tmp2 = zeros(size(tmp));
tmp2(TF) = tmp(TF);
tmp3 = -tmp2;
tmp2(tmp2<0) = 0;
tmp3(tmp3<0) = 0;

[pks1, ~, w1] = findpeaks(tmp2, 'WidthReference', 'halfheight');
[pks2, ~, w2] = findpeaks(tmp3, 'WidthReference', 'halfheight');

ampl = rms([pks1, pks2])/rms(rmoutliers(tmp));
width = mean([w1, w2]);

%% simulate

TRsl = 0.1;
Nsl = 300;
fs=2000;
Ns = TRsl*Nsl*fs;

x_mat = nan(Ns,2,3,3); % samples x ts x corr x spk
x_spks_mat = nan(Ns,2,3,3); % samples x ts x corr x spk
spks_mat = nan(Ns,2,3,3); % samples x ts x corr x spk

ts = ch001(1:58000);
ctr = 1;
spkht = [1 .75 .5];
corr_vals = [.75 .5 .25];
for sh=1:3
    spike_height = spkht(sh);
    for cv=1:3
        corr_val = corr_vals(cv);
        rng(ctr)
        
        % pink noise
        pno = dsp.ColoredNoise('pink',Ns,2); pn = pno();
        
        % white noise
        wno = dsp.ColoredNoise('white',Ns,2); wn = wno();
        
        % alpha noise
        an = bandpass(wno(), [8,12], fs);
        
        % weight noise components
        nf = nan(Ns,2);
        for ii = 1:2
            beta = noisefit(ts, an(:,ii), pn(:,ii), wn(:,ii), fs);
            nf(:,ii) = [an(:,ii) pn(:,ii) wn(:,ii)]*beta; nf(:,ii) = nf(:,ii) / rms(nf(:,ii));
        end

        R = [1 corr_val; corr_val 1];
        L = chol(R);
        nf = nf*L;        
        x = nf;
        x_mat(:,:,cv,sh) = x;
        
        t_ep = (1 : TRsl*fs)';
        locs = [floor(prctile(t_ep, 20)), floor(prctile(t_ep, 80))];
        g = zeros(size(t_ep,1), 1);
        sigma = width/(2*sqrt(2*log(2)));
        for nn = 1:length(locs)
            g_tmp = spike_height*ampl*gaussmf(t_ep, [sigma, locs(nn)]);
            g = g + g_tmp;
            g_apod = g_tmp; g_apod(g_tmp<0.05) = 0;
            base_width = length(g_apod(g_apod~=0));
            g = g - spike_height*ampl*gaussmf(t_ep, [sigma, locs(nn)+base_width]);
        end
        g = repmat(g, Nsl, 1);
        g = g / max(g);
        g = g*5*rms(x);
        
        spks_mat(:,:,cv,sh) = g;
        
        x_spks = x+g;
        x_spks_mat(:,:,cv,sh) = x_spks;

        ctr=ctr+1;
    end
end

x = x_mat;
x_spks = x_spks_mat;
spks = spks_mat;

tmp = (1 : Nsl*TRsl*fs)';
tmp = reshape(tmp, [], Nsl);
tmp2 = tmp(1,:)';
trigs = [tmp2; tmp2(end)+TRsl*fs];

clearvars -except x x_spks spks trigs fs

%%
cd ../..
mkdir tmp

parfor ii = 1 : size(x,2)*size(x,3)*size(x,4)
    [ts,cv,sh] = ind2sub([size(x,2), size(x,3), size(x,4)], ii);
    sf = schrodingerFiltering(x(:,ts,cv,sh), 'fs', fs, 'trigs', trigs, 'passband', [30 700]);
    parsave(sprintf('tmp/%02u_%02u_%02u',ts,cv,sh), sf)
end

sf_mat = nan(size(x));
listing = dir('tmp/*.mat');
for ii = 1:length(listing)
    fname = listing(ii).name;
    load(sprintf('tmp/%s',fname))
    ts = str2double(fname(1:2));
    cv = str2double(fname(4:5));
    sh = str2double(fname(7:8));
    sf_mat(:,ts,cv,sh) = sf;
end
sf = sf_mat;
clear sf_mat
rmdir tmp s

%%
wd = nan(size(x));
at = nan(size(x));
mf = nan(size(x));
for sh = 1:3
    for cv=1:3
        for ts = 1:2
            x_spks_current = x_spks(:,ts,cv,sh);
            dtr = bandpass(x_spks_current, [30 700], fs) ;
            trend = x_spks_current - dtr;
            wd(:,ts,cv,sh) = x_spks_current - wdenoise(dtr);
            mf(:,ts,cv,sh) = medfilt1(dtr) + trend;
            at(:,ts,cv,sh) = filloutliers(dtr,'clip') + trend;
        end
    end
end

%%
corr_err_sf = nan(size(x,3),size(x,4));
corr_err_wd = nan(size(x,3),size(x,4));
corr_err_mf = nan(size(x,3),size(x,4));
corr_err_at = nan(size(x,3),size(x,4));
corr_err_spks = nan(size(x,3),size(x,4));

for sh = 1:3
    for cv=1:3
        corr_gt = corr(x(:,1,cv,sh), x(:,2,cv,sh));

        corr_spks = corr(x_spks(:,1,cv,sh), x_spks(:,2,cv,sh));
        corr_err_spks(cv,sh) = abs(corr_spks - corr_gt);

        corr_sf = corr(sf(:,1,cv,sh), sf(:,2,cv,sh));
        corr_err_sf(cv,sh) = abs(corr_sf - corr_gt);

        corr_wd = corr(wd(:,1,cv,sh), wd(:,2,cv,sh));
        corr_err_wd(cv,sh) = abs(corr_wd - corr_gt);

        corr_at = corr(at(:,1,cv,sh), at(:,2,cv,sh));
        corr_err_at(cv,sh) = abs(corr_at - corr_gt);

        corr_mf = corr(mf(:,1,cv,sh), mf(:,2,cv,sh));
        corr_err_mf(cv,sh) = abs(corr_mf - corr_gt);
    end
end

cd figures/fig10
save('vars.mat')

%%
lbls = categorical({'w. spikes','MF','AT','WD','SF'});
lbls = reordercats(lbls,{'w. spikes','MF','AT','WD','SF'});

for sh = 1:3
    for cv=1:3
        figure
        bar(lbls, [corr_err_spks(cv,sh), corr_err_mf(cv,sh), corr_err_at(cv,sh), corr_err_wd(cv,sh), corr_err_sf(cv,sh)])
        ylabel('correlation error')
        set(gca,'fontsize',18)
    end
end


%%
function beta = noisefit(ts, xa,xp,e, Fs)
% Gabriel Benigno; May 2020
% 
% Weight noise components relative to data given with the EEGSourceSim toolbox (https://doi.org/10.1016/j.jneumeth.2019.108377).
% This code has been slightly modified from the toolbox.
% 

% fprintf('   Weighting noise components');

% load('~/Muller Lab Dropbox/gabriel/schrodinger_filtering/external/EEGSourceSim/+ESSim/+Simulate/private/REC_REO_Averagre_Specs_2sec.mat')

ts = ts(:);
ts = double(ts);
df = Fs/length(ts);
f = (-Fs/2 : df : Fs/2-df)';
FFT = fftshift(abs(fft(ts)));
idx = f<0 | f>100;
f(idx) = [];
FFT(idx) = [];
FFT(f>49 & f<51) = nan;
FFT = fillmissing(FFT,'movmean',2000);

df = Fs/length(xa);
ff = (-Fs/2 : df : Fs/2-df)';
idx = ff<0 | ff>100;
FFTa = fftshift(abs(fft(xa)));
ff(idx) = [];
FFTa(idx) = [];
FFTp = fftshift(abs(fft(xp)));
FFTp(idx) = [];
FFTe = fftshift(abs(fft(e)));
FFTe(idx) = [];

FFT = interp1(f, FFT, ff);
FFT(isnan(FFT)) = 0;

beta = [FFTa FFTp FFTe] \ FFT;

% plot(ff,FFT, ff, [FFTa FFTp FFTe]*beta )

% fprintf('...Done \n')

end


function parsave(fname, sf)
save(fname, 'sf')
end
