clearvars;clc
load('../../realdata/data1/trigs.mat')
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
cd ../..
mkdir tmp
parfor ch = 1:30
    sf = schrodingerFiltering(aas{ch}, 'fs', 2048, 'mainsfreq',50, 'trigs',trigs,'passband',[30, 700]);
    parsave(sprintf('tmp/%02u.mat',ch), sf)
end
%%
listing = dir('tmp/*.mat');
sf_cell = cell(30,1);
for ii = 1:length(listing)
    fname=sprintf('tmp/%s',listing(ii).name);
    load(fname)
    ch = str2double(fname(5:6));
    sf_cell{ch} = sf;
end
sf = sf_cell;
clear sf_cell
rmdir tmp s

%% just the 7 evoked channels
figure
idxs = [14 15 16 18 19 26 27];
os=0;
dos=300;
tt=trigs(1) : trigs(end)-1;
for ii = 1:7
    plot(tt/2048,aas{idxs(ii)}(tt)+os,'k',   tt/2048, sf{idxs(ii)}(tt)+os,'b','linewidth',0.1)
    os=os+dos;
    hold on
end
hold off
xlim([tt(1),tt(end)]/2048)
ylim([min(aas{1}(tt)), max(aas{27}(tt)+6*dos)])
set(gca,'ytick',[])
set(gca,'ticklength',[0 0])
set(gca,'fontsize',26)
xlabel('time (s)','FontSize',32)


%% CWT
nmbr=5;
[wt,f]=cwt(aas{idxs(nmbr)}(tt),2048);

figure(1)
imagesc(tt/2048,f,  abs(wt))
set(gca,'YDir','normal','yscale','log','fontsize',18)
xlabel('time (s)')
ylabel('frequency (Hz)')
ylim([1 100])
caxis([0 100])

[wt,f]=cwt(sf{idxs(nmbr)}(tt),2048);

figure(2)
imagesc(tt/2048,f,  abs(wt))
set(gca,'YDir','normal','yscale','log','fontsize',18)
xlabel('time (s)')
ylabel('frequency (Hz)')
ylim([1 100])
caxis([0 100])

%% timeseries closeup and psd
ch=1;
tt=trigs(1) : trigs(end)-1;
[pxx,f]=pwelch(aas{ch}(tt),[],[],[],2048);
pxx2=pwelch(sf{ch}(tt),[],[],[],2048);
ref1 = pxx; ref2 = pxx2;
for ii = 1:100
    ref1 = ref1(randperm(length(f)));
    ref2 = ref2(randperm(length(f)));
end

P = 10*log10(pxx./ref1);
P2 = 10*log10(pxx2./ref2);

figure(1)
semilogx(f,P,'k',f,P2,'r','linewidth',2)
xlim([f(1),100])
xlabel('frequency (Hz)')
ylabel('power (dB)')
legend('AAS','AAS+SF')

figure(2)
plot(tt/2048,aas{ch}(tt),'k',tt/2048,sf{ch}(tt),'r')
ttt = find(tt>=2048*33 & tt<=2048*36);
xlabel('time (s)')
ylabel('voltage (µV)')
legend('AAS','AAS+SF')





%%
function parsave(fname, sf)
save(fname, 'sf')
end
