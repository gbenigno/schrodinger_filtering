clearvars;clc
load('../../realdata/data1/raw/ch007.mat')
raw = ch007;
load('../../realdata/data1/aas/ch007.mat')
aas = ch007;
load('../../realdata/data1/trigs.mat')
fs = 2048;
load('../../realdata/data1/times.mat')

xlims=[trigs(1), trigs(7)]/fs;
width=1000; 
height = 200;
linewidthh=0.5;
fontsizee=12; 

cutoff = 400;

%% 1A: six epochs of raw
figure
plot(times,raw,'k','linewidth',linewidthh)
ylabel('voltage (µV)')
xlim(xlims)
set(gca,'fontsize',fontsizee)
set(gca,'ticklength',[0 0])
pbaspect([5 1 1])

%% 1B: closeup of six raw epochs; AAS; AAS-->LPF
figure
plot(times, raw, 'k', ...
    times, lowpass(aas, cutoff, fs), 'r','LineWidth',1)
ylims=[-130 130];
xlabel('time (s)')
ylabel('voltage (µV)')
xlim(xlims)
ylim(ylims)
set(gca,'fontsize',fontsizee)
set(gca,'ticklength',[0 0])
pbaspect([5 1 1])

%% 1C: PSDs of raw, AAS
figure
tt=trigs(1):trigs(end)-1;
Ns = length(tt);
y = raw(tt);
yaas = lowpass(aas(tt), cutoff, fs);

Y = pwelch(y, [], [], [], fs);
[Yaas, f] = pwelch(yaas, [], [], [], fs);

semilogy(f,Y,'k',f,Yaas,'r');

xl1 = 0; xl2 = 160;
xlims = [xl1 xl2];

yl1 = min([Y(f>=xl1 & f<=xl2); Yaas(f>=xl1 & f<=xl2)]);
yl2 = max([Y(f>=xl1 & f<=xl2); Yaas(f>=xl1 & f<=xl2)])+10;

ylims=[yl1,yl2];
xlim(xlims)
ylim(ylims)
set(gca,'ticklength',[0 0])
xticks([0 7:7:7*22])
xticklabels({'0','7','14','21','28','35','42','49','56','63','70','77','84','91','98','105','112','119','126','133','140','147','154'})
xlabel('frequency (Hz)')
ylabel('power (µV^2/Hz)')
pbaspect([5 1 1])
set(gca,'fontsize',10)