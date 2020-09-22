clearvars;clc

% get given positions of electrodes
% Get x, y, and z positions of electrodes and cnonvert to spherical coordinates
load('../../realdata/data1/chanlocs.mat')

x = [chanlocs.X];
y = [chanlocs.Y];
z = [chanlocs.Z];
v = [x; y; z];
v0 = [(x(18)+x(19))/2; (y(18)+y(19))/2; (z(18)+z(19))/2];
d=acos(sum(v0 .* v));
af = exp(-d.^2);
af = af/max(af)*2;

TRsl = 0.1;
Nsl = 600;
fs=2000;

t = (1 : Nsl*TRsl*fs)';
Ns = length(t);

% pink noise
pno = dsp.ColoredNoise('pink',Ns,1); pn = pno();

% white noise
wno = dsp.ColoredNoise('white',Ns,1); wn = wno();
wno_spurt = dsp.ColoredNoise('white', 5*fs, 1);

gamma_act = zeros(size(t));
for kk = [5 15 25 35 45 55]*fs
    gamma_it = zeros(size(t));
    gamma_it(kk) = 1;
    gamma_act = gamma_act + conv(gamma_it, bandpass(wno_spurt(), [40,60], fs), 'same');
end
gamma_act = gamma_act/rms(gamma_act(gamma_act~=0));
gamma_act = af(15)*gamma_act;

% alpha noise
an = bandpass(wno(), [8,12], fs);

% weight noise components
load('../../realdata/data1/raw/ch015.mat')
ts = ch015(1:58000);
beta = noisefit(ts, an, pn, wn, fs);
nf = [an pn wn]*beta; nf = nf / rms(nf); 

x = nf + gamma_act;

%%
figure
plot(t,pn,t,wn,t,an,t,gamma_act)
lgd=legend('pink noise','white noise','alpha noise','gamma activity','Orientation','horizontal');
xlim([0 6e4])
set(gca,'ticklength',[0 0])
set(gca,'ytick',[])
set(gca,'fontsize',16)
lgd.FontSize = 12;
xlabel('time (s)')

%%
figure
eeglabpath = '../../external/eeglab2019_1';
addpath(genpath(eeglabpath))
topoplot(af,chanlocs,'whitebk','on');
c=colorbar;
c.Label.String = 'gamma activity';
rmpath(genpath(eeglabpath))






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


