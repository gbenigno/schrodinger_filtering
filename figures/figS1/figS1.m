clearvars;clc
%%
ch=19;
names = {'ground truth','SF','WD','AT','MF'};
spkampl = {'high','medium','low'};
set(gcf,'color','white')
for ii = 1:3 % spike height
    for jj = 1:5 % technique
        figure
        switch jj
            case 1
                cwt(x(:,ch,ii),fs)
            case 2
                cwt(sf(:,ch,ii),fs)
            case 3
                cwt(wd(:,ch,ii),fs)
            case 4
                cwt(at(:,ch,ii),fs)
            case 5
                cwt(mf(:,ch,ii),fs)
        end
        title(sprintf('%s; %s spike amplitudes',names{jj},spkampl{ii}))
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end