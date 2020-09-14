TRsl = 0.1; pausepercent = 0.2; stf = 700; nch = 1; Nsl = 0; fs=2048;
addpath '~/Muller Lab Dropbox/gabriel/schrodinger_filtering'
addpath '~/Muller Lab Dropbox/gabriel/schrodinger_filtering/simulation'
addpath '~/Muller Lab Dropbox/gabriel/schrodinger_filtering/simulation/spike_simulation/'

sar_sf = nan(1000,1);
sar_at = nan(1000,1);
sar_mf = nan(1000,1);
sar_wd = nan(1000,1);

sf_record = nan(408,1000);
x_record = nan(408,1000);
spks_record = nan(408,1000);

[t, ~, ~, ~, Ns] = gradsim(fs, TRsl, Nsl, pausepercent, stf);

%%
parfor ii = 1:1000
    
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
    
end
%%

sf=sf_record;
x = x_record;
spks = spks_record;
clear sf_record x_record spks_record

save('~/Muller Lab Dropbox/gabriel/schrodinger_filtering/figures/fig6/wkspc.mat')

%

%{
plot(sar_wd,sar_sf,'o',sar_at,sar_sf,'o',sar_mf, sar_sf, 'o', [0, max([sar_wd; sar_at; sar_mf])], [0, max([sar_wd; sar_at; sar_mf])], [0, max([sar_wd; sar_at; sar_mf])], 4*[0, max([sar_wd; sar_at; sar_mf])])
ylabel('signal-to-artifact ratio of SF');xlabel('signal-to-artifact ratio of other methods')
% legend('comparison to WD','comparison to AT','y=x')
pbaspect([3 1 1])
xlim([0 max(sar_wd)*1.1]);ylim([0 max(sar_sf)*1.1]);
% title('signal-to-artifact ratio')
%}

%asfsad