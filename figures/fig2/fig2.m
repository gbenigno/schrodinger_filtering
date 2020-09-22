clearvars;clc
TRsl = 0.1; pausepercent = 0.2; stf = 700; nch = 1; Nsl = 100;
fs = 2048; snrr = 1; sh = 4;
[t, ~, trigs_orig, ~, Ns] = gradsim(fs, TRsl, Nsl, pausepercent, stf); % gradient artifact
trigs = [trigs_orig; trigs_orig(end) + min(diff(trigs_orig))];
pno = dsp.ColoredNoise('pink',Ns,nch);
wno = dsp.ColoredNoise('white',Ns,nch);
s = ssvep(t);
pn = pno();
wn = wno();
an = bandpass(wno(), [8,12], fs);

load('../../realdata/data1/raw/ch001.mat')
ts = ch001(1:58000);
beta = noisefit(ts, an, pn, wn, fs);
nf = [an pn wn]*beta; nf = nf / rms(nf);
x =  signoise(snrr, s, nf);

g1 = 2*rms(x)*gaussmf(t,[8/2048, t(randi([1, length(t)]))]);
g2 = 2*rms(x)*gaussmf(t,[8/2048, t(randi([1, length(t)]))]);
g3 = -2*rms(x)*gaussmf(t,[8/2048, t(randi([1, length(t)]))]);

f = -fs/2: 1/t(end) : fs/2;

G1 = abs(fftshift(fft(g1)));
G2 = abs(fftshift(fft(g2)));
G3 = abs(fftshift(fft(g3)));

figure
tl=tiledlayout(2,3);
nexttile(1)
plot(t,g1,'k')
xlim([t(1),t(end)])
ylim([-1,4])
set(gca,'ticklength',[0,0])
xlabel('time (s)')

nexttile(2)
plot(t,g2,'k')
xlim([t(1),t(end)])
ylim([-1,4])
set(gca,'ticklength',[0,0])
xlabel('time (s)')

nexttile(3)
plot(t,g3,'k')
xlabel('time (s)')
tl.TileSpacing='compact';
xlim([t(1),t(end)])
ylim([-4,1])
set(gca,'ticklength',[0,0])

nexttile([1 3])
plot(f,G1,'k')
set(gca,'ticklength',[0,0])
xlim([f(1),f(end)])
xlabel('frequency (Hz)')

%% bunch of random spikes look like pink noise

g = zeros(size(t));
for ii = 1:20
    r = round(rand); if r ==0, r=-1; end
    g = g + r*6*rms(x)*gaussmf(t,[8/2048, t(randi([1, length(t)]))]);
end

G = abs(fftshift(fft(g)));
X = abs(fftshift(fft(x)));

figure
tl=tiledlayout(2,1);

nexttile
plot(t,x+g,'k')
xlim([t(1),t(end)])
set(gca,'ticklength',[0,0])
xlabel('time (s)')
legend('EEG+spikes')

nexttile
plot(f,X,f,G)
set(gca,'ticklength',[0,0])
xlim([-120,120])
ylim([0,max(G)])
xlabel('frequency (Hz)')
legend('EEG','spikes')
tl.TileSpacing='compact';

%% --------------------------------------------------------
function s = ssvep(t)
a1 = 2*rand;
a2 = 1.5*rand;
f1 = 2*rand;
f2 = 6*rand;
ph1 = pi*rand - pi/2;
ph2 = pi*rand - pi/2;
s1 = a1*cos(2*pi*f1*t - ph1);
s2 = a2*cos(2*pi*f2*t - ph2);
s = s1+s2;
end


function [x,s] =  signoise(snrr, s, nf)
m = snrr * (rms(nf)/rms(s))^2;
s = m*s;
x = s+nf;
end

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
