clearvars;clc

load('../../realdata/data1/chanlocs.mat')
aas = cell(30,1);
raw = cell(30,1);
for ch = 1:30
    load(sprintf('../../realdata/data1/aas/ch%03u.mat',ch))
    eval(sprintf('aas{ch} = ch%03u;',ch))
    load(sprintf('../../realdata/data1/raw/ch%03u.mat',ch))
    eval(sprintf('raw{ch} = ch%03u;',ch))
    clear(sprintf('ch%03u',ch))
end

%%
% get given positions of electrodes
% Get x, y, and z positions of electrodes and cnonvert to spherical coordinates
x = [chanlocs.X];
y = [chanlocs.Y];
z = [chanlocs.Z];
v = [x; y; z];
v0 = [(x(18)+x(19))/2; (y(18)+y(19))/2; (z(18)+z(19))/2];
d=acos(sum(v0 .* v));
af = exp(-d.^2);
af = af/max(af)*2;

TRsl = 0.1;
Nsl = 300;
fs=2000;
Ns = TRsl*Nsl*fs;

x_mat = nan(Ns,30,3);
spks_mat = nan(Ns,30,3);
x_spks_mat = nan(Ns,30,3);
alpha_act_mat = nan(Ns,30,3);
trend_mat = nan(Ns,30,3);

tt = 58887 : 304369+300;
tmp = double(aas{7}(tt));
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
parfor ch = 1:30
    
    fctrs = [1 .75 .5];
    t = (1 : Ns)';
    
    for ii = 1:3
        rng(sub2ind([30,3],ch,ii))
        
        fctr = fctrs(ii);
        
        % pink noise
        pno = dsp.ColoredNoise('pink',Ns,1); pn = pno();
        
        % white noise
        wno = dsp.ColoredNoise('white',Ns,1); wn = wno();
        wno_spurt = dsp.ColoredNoise('white', 5*fs, 1);
        
        alpha_act = zeros(size(t));
        for kk = [5 15 25]*fs
            alpha_it = zeros(size(t));
            alpha_it(kk) = 1;
            alpha_act = alpha_act + conv(alpha_it, bandpass(wno_spurt(), [8,12], fs), 'same');
        end
        alpha_act = alpha_act/rms(alpha_act(alpha_act~=0));
        alpha_act = af(ch)*alpha_act;
        
        % alpha noise
        an = bandpass(wno(), [8,12], fs);

        % weight noise components
        ts = raw{ch}(1:58000);
        beta = noisefit(ts, an, pn, wn, fs);
        nf = [an pn wn]*beta; nf = nf / rms(nf); 
        
        x = nf + alpha_act;
        alpha_act_mat(:,ch,ii) = alpha_act;
        x_mat(:,ch,ii) = x;
        
        t_ep = (1 : TRsl*fs)';
        locs = [floor(prctile(t_ep, 20)), floor(prctile(t_ep, 80))];
        g = zeros(size(t_ep,1), 1);
        sigma = width/(2*sqrt(2*log(2)));
        for nn = 1:length(locs)
            g_tmp = fctr*ampl*gaussmf(t_ep, [sigma, locs(nn)]);
            g = g + g_tmp;
            g_apod = g_tmp; g_apod(g_tmp<0.05) = 0;
            base_width = length(g_apod(g_apod~=0));
            g = g - fctr*ampl*gaussmf(t_ep, [sigma, locs(nn)+base_width]);
        end
        g = repmat(g, Nsl, 1);

        spks_mat(:,ch,ii) = g;
        x_spks = x+g;   
        x_spks_mat(:,ch,ii) = x_spks;
        
        trend = x_spks - bandpass(x_spks, [30 700], fs);
        trend_mat(:,ch,ii) = trend;
        
    end
end

x = x_mat;
spks = spks_mat;
x_spks = x_spks_mat;
alpha_act = alpha_act_mat;
trend = trend_mat;

d1 = TRsl*fs;
x_rs = reshape(x, d1, [], 30, 3);
spks_rs = reshape(spks, d1, [], 30, 3);
x_spks_rs = reshape(x_spks, d1, [], 30, 3);
alpha_act_rs = reshape(alpha_act, d1, [], 30, 3);
trend_rs = reshape(trend, d1, [], 30, 3);

tmp = (1 : Nsl*TRsl*fs)';
tmp = reshape(tmp, [], Nsl);
tmp2 = tmp(1,:)';
trigs = [tmp2; tmp2(end)+TRsl*fs];

%%
cd ../..
mkdir tmp
parfor ii = 1:90
    [ch,jj] = ind2sub([30,3],ii);
    fprintf('%u\n',ii)
    sf = schrodingerFiltering(x_spks(:,ch,jj), 'fs', fs, 'passband', [30 700], 'trigs', trigs);
    parsave(sprintf('tmp/%02u_%02u',ch,jj), sf)
end

%%
listing = dir('tmp/*.mat');
sf_mat = nan(size(x));
for jj = 1:length(listing)
    fname = sprintf('tmp/%s',listing(jj).name);
    load(fname)
    ch = str2double(fname(5:6));
    ii = str2double(fname(8:9));
    sf_mat(:,ch,ii) = sf;
end
sf = sf_mat;
sf_rs = reshape(sf, [], Nsl, size(x,2), size(x,3));

cd figures/fig8
save('vars.mat')

rmdir ../../tmp s

%%
wd = nan(size(x));
at = nan(size(x));
mf = nan(size(x));

sar_sf = nan(size(x));
sar_wd = nan(size(x));
sar_mf = nan(size(x));
sar_at = nan(size(x));

for ch=1:30
    for ii=1:3
        dtr = x_spks(:,ch,ii) - trend(:,ch,ii);
        wd(:,ch,ii) = x_spks(:,ch,ii) - wdenoise(dtr);
        mf(:,ch,ii) = medfilt1(dtr) + trend(:,ch,ii);
        at(:,ch,ii) = filloutliers(dtr,'clip') + trend(:,ch,ii);
    end
end

wd_rs = reshape(wd, size(x_rs,1), [], size(x_rs,3), size(x_rs,4));
at_rs = reshape(at, size(x_rs,1), [], size(x_rs,3), size(x_rs,4));
mf_rs = reshape(mf, size(x_rs,1), [], size(x_rs,3), size(x_rs,4));

for ch = 1:30
    for ii = 1:3
        for ep = 1:Nsl
            sar_sf(ep,ch,ii) = sum( (sf_rs(:,ep,ch,ii) - spks_rs(:,ep,ch,ii)).^2) ./ sum((sf_rs(:,ep,ch,ii) - x_rs(:,ep,ch,ii)).^2);
            sar_wd(ep,ch,ii) = sum((wd_rs(:,ep,ch,ii) - spks_rs(:,ep,ch,ii)).^2) ./ sum((wd_rs(:,ep,ch,ii) - x_rs(:,ep,ch,ii)).^2);
            sar_at(ep,ch,ii) = sum((at_rs(:,ep,ch,ii) - spks_rs(:,ep,ch,ii)).^2) ./ sum((at_rs(:,ep,ch,ii) - x_rs(:,ep,ch,ii)).^2);
            sar_mf(ep,ch,ii) = sum((mf_rs(:,ep,ch,ii) - spks_rs(:,ep,ch,ii)).^2) ./ sum((mf_rs(:,ep,ch,ii) - x_rs(:,ep,ch,ii)).^2);
        end
    end
end

%% high
sar_wd_tmp = sar_wd(:,:,1); sar_wd_tmp = sar_wd_tmp(:);
sar_mf_tmp = sar_mf(:,:,1); sar_mf_tmp = sar_mf_tmp(:);
sar_sf_tmp = sar_sf(:,:,1); sar_sf_tmp = sar_sf_tmp(:);
sar_at_tmp = sar_at(:,:,1); sar_at_tmp = sar_at_tmp(:);

scatter(sar_wd_tmp, sar_sf_tmp)
hold on
scatter(sar_at_tmp, sar_sf_tmp)
scatter(sar_mf_tmp, sar_sf_tmp)
plot(1:100, 1:100, '--k',  1:100, 4*(1:100), '--k','LineWidth',3)

p1 = plot(nan,nan,'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerSize',10);
p2 = plot(nan,nan,'o', 'color', [0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',10);
p3 = plot(nan,nan,'o', 'color', [0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerSize',10);

hold off
set(gca,'yscale','log','xscale','log','fontsize',28)
legend([p1 p2 p3],'comparison to WD','comparison to AT','comparison to MF','Location','southeast')
pbaspect([1 1 1])
xlabel('signal-to-artifact ratio of other technique')
ylabel('signal-to-artifact ratio of SF')

%% med

set(gcf,'color','white')

sar_wd_tmp = sar_wd(:,:,2); sar_wd_tmp = sar_wd_tmp(:);
sar_mf_tmp = sar_mf(:,:,2); sar_mf_tmp = sar_mf_tmp(:);
sar_sf_tmp = sar_sf(:,:,2); sar_sf_tmp = sar_sf_tmp(:);
sar_at_tmp = sar_at(:,:,2); sar_at_tmp = sar_at_tmp(:);

scatter(sar_wd_tmp, sar_sf_tmp)
hold on
scatter(sar_at_tmp, sar_sf_tmp)
scatter(sar_mf_tmp, sar_sf_tmp)
plot(1:100, 1:100, '--k',1:100, 4*(1:100), '--k', 'LineWidth',3)

p1 = plot(nan,nan,'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerSize',10);
p2 = plot(nan,nan,'o', 'color', [0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',10);
p3 = plot(nan,nan,'o', 'color', [0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerSize',10);

hold off
set(gca,'yscale','log','xscale','log','fontsize',28)
legend([p1 p2 p3],'comparison to WD','comparison to AT','comparison to MF','Location','southeast')
pbaspect([1 1 1])
xlabel('signal-to-artifact ratio of other technique')
ylabel('signal-to-artifact ratio of SF')


%% low

set(gcf,'color','white')

sar_wd_tmp = sar_wd(:,:,3); sar_wd_tmp = sar_wd_tmp(:);
sar_mf_tmp = sar_mf(:,:,3); sar_mf_tmp = sar_mf_tmp(:);
sar_sf_tmp = sar_sf(:,:,3); sar_sf_tmp = sar_sf_tmp(:);
sar_at_tmp = sar_at(:,:,3); sar_at_tmp = sar_at_tmp(:);

scatter(sar_wd_tmp, sar_sf_tmp)
hold on
scatter(sar_at_tmp, sar_sf_tmp)
scatter(sar_mf_tmp, sar_sf_tmp)
plot(1:100, 1:100, '--k',1:100, 4*(1:100), '--k',  'LineWidth', 3)

p1 = plot(nan,nan,'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerSize',10);
p2 = plot(nan,nan,'o', 'color', [0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',10);
p3 = plot(nan,nan,'o', 'color', [0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerSize',10);

hold off
set(gca,'yscale','log','xscale','log','fontsize',28)
legend([p1 p2 p3],'comparison to WD','comparison to AT','comparison to MF','Location','southeast')
pbaspect([1 1 1])
xlabel('signal-to-artifact ratio of other technique')
ylabel('signal-to-artifact ratio of SF')



%% topoplots
t = (1:18e3)'/fs;
idx = (t>=2.5 & t<=7.5) | (t>=12.5 & t<=17.5) | (t>=22.5 & t<=27.5) | (t>=32.5 & t<=37.5) | (t>=42.5 & t<=47.5) | (t>=52.5 & t<=57.5);

sh = {'high','med','low'};
tc = {'gt','sf'};
titles = {'Ground Truth','SF'};

close all

eeglabpath = '../..external/eeglab2019_1';
addpath(genpath(eeglabpath))

aa = nan(30,3);
for jj = 1:3
    for kk = 1:2
        for ii = 1:30
            switch kk
                case 1
                    aa(ii,jj) = rms(bandpass(x(:,ii,jj),[8 12],fs));
                case 2
                    aa(ii,jj) = rms(bandpass(sf(:,ii,jj),[8 12],fs));
            end
        end
        
        topoplot(aa(:,jj),EEG.chanlocs,'whitebk','on');
        if kk==1
            c=colorbar; c.Label.String = 'alpha activity'; c.FontSize = 18;
        end
        title(titles{kk},'FontSize',24)
        keyboard
        close all
    end
end

rmpath(genpath(eeglabpath))

%% find exemplary epochs
names={'ts_high.pdf','ts_med.pdf','ts_low.pdf'};
t = (1:200)'/fs;
for ii = 1:3
    ratio = sar_sf(:,:,ii) ./ sar_wd(:,:,ii);
    ratio = ratio(:);
    [~,I] = max(ratio);
    [r,c] = ind2sub([600 30], I);
    
    for jj = 1:5
        switch jj
            case {1}
                plot(t, x_spks_rs(:,r,c,ii) + 10*(jj-1)  ,'k', t, x_rs(:,r,c,ii)+10*(jj-1), 'r')
                hold on
            case {2}
                plot(t, x_spks_rs(:,r,c,ii) + 10*(jj-1)  ,'k', t, sf_rs(:,r,c,ii)+10*(jj-1), 'r')
            case {3}
                plot(t, x_spks_rs(:,r,c,ii) + 10*(jj-1)  ,'k', t, wd_rs(:,r,c,ii)+10*(jj-1), 'r')
            case {4}
                plot(t, x_spks_rs(:,r,c,ii) + 10*(jj-1)  ,'k', t, at_rs(:,r,c,ii)+10*(jj-1), 'r')
            case {5}
                plot(t, x_spks_rs(:,r,c,ii) + 10*(jj-1)  ,'k', t, mf_rs(:,r,c,ii)+10*(jj-1), 'r')
        end  
    end
    ylim([-7 47])
    xlabel('time (s)')
    pbaspect([1 1 1])
    hold off
    set(gcf,'color','none')
    export_fig(names{ii})
    
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