clearvars;clc

%% AAS
%{
for ii = 1:17
    fprintf('ii=%u, ',ii)
    cd(sprintf('/Volumes/GB2TB/osfdataset/sub-%03u/sourcedata-eeg_inside-MRT/eeg',ii))
    listing = dir('*.mat');
    for jj = 1:length(listing)
        fprintf('jj=%u\n',jj)
        [~,fname,~] = fileparts(listing(jj).name);
        load(listing(jj).name)
        EEG_aas = AAS(EEG);
        save(sprintf('/Volumes/GB2TB/osfdataset/gb/inside/aas/sub-%03u/%s',ii,fname),'EEG_aas')
    end
end
%}

%% timeseries
%{
load('/Volumes/GB2TB/osfdataset/gb/inside/aas/sub-017/sub-017_task-pdm_acq-insideMRT_run-03_eeg.mat')
t = EEG_aas.times/1000;

figure
plot(t,EEG_aas.data(1,:))
ylim([-600 600])
xlim([t(1) t(end)])
pbaspect([2 1 1])
xlabel('time (s)')
ylabel('voltage (µV)')

ind = t>481.6 & t<482.6;
tt = t(ind);
y = EEG_aas.data(1,ind);
trigs = (1:1000:length(y)+1)';
cd ../..
ysf = schrodingerFiltering(y,'fs',5000,'passband',[30, 700],'trigs',trigs);

figure
plot(tt,y,tt,ysf)
xlabel('time (s)')
ylabel('voltage (µV)')
%}

%% epoching, despiking epochs, and averaging
%{
clearvars; clc
addpath ../..
fs = 5000;
cd('/Volumes/GB2TB/osfdataset/gb/inside/aas')
listing = dir('sub*');
erp = [];%cell(4,1); % four different conditions
erp_ch = cell(64,1);
cdns = {'S 10', 'S 20', 'S 30', 'S 40'; 'S 11', 'S 21', 'S 31', 'S 41'}'; % left/right for each condition
for fldr = 1:length(listing)
    fprintf('%u, ',fldr)
    subj = str2double(listing(fldr).name(end-2 : end));
    switch subj
        case {1}
            runs = [1 2 3 4];
        case {3}
            runs = [1 3 4];                                        
        case {6}
            runs = [1 2 3 5];  
        case {13, 14, 15, 16, 17}
            runs = [1 3 4 5]; 
        case {5, 7, 8, 9, 11, 12}
            runs = [1 2 3 4 5]; 
        otherwise
            runs = [1 2 3 4 5]; 
    end
    for run = runs
        fname = sprintf('sub-%03u/sub-%03d_task-pdm_acq-insideMRT_run-%02u_eeg.mat', subj, subj, run);
        load(fname)
        trigs = gettrigs(EEG_aas);
        data = EEG_aas.data(:,trigs(1):trigs(end));
        avg = mean(data);
        data_rr = data - avg;
        for ii = [10 46 60 9 45 59]
            data_rr(ii,:) = fftgausshp(data_rr(ii,:),1,fs);
        end
        for LR = 1:2
            switch LR
                case {1}
                    chs = [10; 46; 60];
                case {2}
                    chs = [9; 45; 59];
            end
            for cdn = 3%1:size(cdns,1)
                idxs = epoching(EEG_aas.event, cdns{cdn, LR}) - trigs(1) + 1;
                for jj = 1:length(idxs)
                    idx = idxs(jj);
                    for ii = 1:length(chs)
                        ch = chs(ii);
                        try
                            tmp = data_rr(ch, idx-100e-3*fs+1 : idx+500e-3*fs )';
                            erp = [erp, tmp];
                        catch
                        end
                    end
                    
                    switch subj
                        case {1,2,3,4,5}
                            for pp = [1:30, 33:64]
                                try
                                    tmp = data_rr(pp, idx-100e-3*fs+1 : idx+500e-3*fs )';
                                    erp_ch{pp} = [erp_ch{pp}, tmp];
                                catch
                                end
                            end
                    end
                    
                end
            end
        end
    end
end
erp_avg = nanmean(erp, 2);
erp_ch_avg = nan(3000,64);
for ii = [1:30, 33:64]
    erp_ch_avg(:,ii) = nanmean(erp_ch{ii},2);
end
cd('~/Muller Lab Dropbox/gabriel/schrodinger_filtering/figures/fig9')
save('erp.mat','erp','erp_avg','erp_ch','erp_ch_avg')
%}

%% SF overall
%{
clearvars;clc
load('erp.mat')

addpath ../..

erp_sf = nan(3000,4104);

parfor ii = 1 : 4104
    try
        erp_sf(:,ii) = schrodingerFiltering(double(erp(:,ii)),'fs',5000,'trigs',(1:100:3001)','passband',[30 700] );
%         dummy='';
%         parsave(sprintf('sf/%u.txt',ii),dummy)
    catch
    end
end
erp_sf_avg = nanmean(erp_sf,2);
save('erp.mat','erp_sf','erp_sf_avg','-append')

fs=5000;
ind = -100e-3*fs+1 : 500e-3*fs;

figure
plot(ind/fs, erp_avg, ind/fs, erp_sf_avg)
xlabel('time relative to stimulus (s)')
ylabel('voltage (µV)')
%}

%% SF per channel

clearvars;clc
load('erp.mat')
tmp = erp_ch';
tmp{31} = single(nan(3000,size(tmp{1},2)));
tmp{32} = single(nan(3000,size(tmp{1},2)));
tmp = cell2mat(tmp);

erp_ch_sf = nan(size(tmp));

addpath ../..

parfor ii = 1:size(tmp,2)
    try
        erp_ch_sf(:,ii) = schrodingerFiltering(double(tmp(:,ii)),'fs',5000,'trigs',(1:100:3001)','passband',[30 700] );
%         dummy='';
%         parsave(sprintf('sf/%u.txt',ii),dummy)
    catch
    end
end
erp_ch_sf = mat2cell(erp_ch_sf, 3000, repelem(size(erp_ch_sf,2)/64,1,64))';
erp_ch_sf{31} = [];
erp_ch_sf{32} = [];
erp_ch_sf_avg = nan(3000,64);
for ii = [1:30, 33:64]
    erp_ch_sf_avg(:,ii) = nanmean(erp_ch_sf{ii},2);
end
save('erp.mat','erp_ch_sf','erp_ch_sf_avg','-append')
%}


%%

clearvars;clc
load('chanlocs.mat')
load('erp.mat','erp_ch_avg','erp_ch_sf_avg')
fs=5000;
time = (1/fs : 1/fs : 3000/fs)' - 0.1;

eeglabpath='../../external/eeglab2019_1';
addpath(genpath(eeglabpath))
frm = 2300;
x = erp_ch_avg(frm, :);
xsf = erp_ch_sf_avg(frm, :);

figure

topoplot(x, chanlocs,  'maplimits', [-5, 5])%,'shading','interp')
c = colorbar;
c.Label.String = 'voltage (µV)';
c.FontSize=18;
title('AAS','FontSize',20)

figure
topoplot(xsf, chanlocs,  'maplimits', [-5, 5])%,'shading','interp')
c = colorbar;
c.Label.String = 'voltage (µV)';
c.FontSize=18;
title('AAS+SF','FontSize',20)

rmpath(genpath(eeglabpath))
%}

%%
function parsave(fname, dummy)
save(fname, 'dummy')
end

function trigs = gettrigs(EEGin)

trigs = [];
ctr=0;
for ii = 1:length(EEGin.event)
    if strcmp(EEGin.event(ii).type, 'R128')
        ctr=ctr+1;
        if isempty(trigs)
            trigs = [trigs; EEGin.event(ii).latency];
        end
    end
end
for ii = 1:ctr
    trigs = [trigs; trigs(end)+10000];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTGAUSSHP  Gauss Filter with FFT
%
%  Copyright (c) 2010-2012 Johann Glaser <Johann.Glaser@gmx.at>
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
% MA  02110-1301, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fftgausshp(in,fc,srate)
  % FFTGAUSSHP  Gauss Filter with FFT
  %
  % Implements a filter in frequency domain by multiplying the input
  % signal's FFT with the transfer function.
  %
  %  in       signal to be filtered
  %  fc       3dB cut-off frequency relative to nyquist frequency
  %  out      result
  %
  % This function fully removes the DC part, so no extra DC removal is
  % required.

  % Perform FFT
  X=fft(in);

  % gauss curve parameter
  % H = 1-exp(-(f/2s)^2)
  % f = ?
  % 1-H = exp(-(f/2s)^2)
  % ln(1-H) = -(f/2s)^2
  % sqrt(-ln(1-H)) = f/2s
  % f = 2s*sqrt(-ln(1-H))
  %
  % H = 1/sqrt(2) -> f = ?
  % f = 2 * s * 1.1081
  % f = 2.2163 * s;
  %
  % s = ?
  % s = f / (2*sqrt(-ln(1-H)))
  % s = f / 2.2163
  % s = f * 0.4512
  sigma = fc / (2*sqrt(-log(1-1/sqrt(2))));

  % transfer function
  nyq = floor(length(X)/2);
  sigma = sigma * nyq/(srate/2);
  f = 0:nyq;
  H = 1.0 - exp(-(f / (2 * sigma)).^2);
  % Note: H(1) = 0.0, i.e. full DC block

  % apply transfer function
  % f(1) ..... 0Hz
  % f(2) ..... srate/length(in) Hz
  % ...
  % f(end) ... -srate/length(in) Hz
  if size(in,1) == 1   % single-row vector
    X(1:nyq)           = X(1:nyq)           .* H(1:nyq);
    X((end-nyq+2):end) = X((end-nyq+2):end) .* H(nyq:-1:2);
    if rem(length(X),2)==0
      % even signal length: set unreached point to H(nyq)
      X(nyq+1) = H(nyq);
    end
  else
    % fft() is performed to each column (i.e. 2nd index) along the row
    % (i.e. 1st index). Therefore we have to set
    for i=1:size(X,2)
      X(1:nyq,i)           = X(1:nyq,i)           .* H(1:nyq)';
      X((end-nyq+2):end,i) = X((end-nyq+2):end,i) .* H(nyq:-1:2)';
      if rem(length(X),2)==0
        % even signal length: set unreached point to H(nyq)
        X(nyq+1,i) = H(nyq);
      end
    end
  end

  % Inverse FFT
  out = ifft(X);

end


function M = epoching(strct, typee) 
M = [];
for ii = 1:length(strct)
    if strcmp(strct(ii).type, typee)
        M = [M; strct(ii).latency];
    end
end
end


function EEGout = AAS(EEGin)

nch = EEGin.nbchan;
fs = EEGin.srate;

trigs = gettrigs(EEGin);
avgmat = getavgmat(trigs);
EEGout = EEGin;
EEGin_data = EEGin.data;
EEGout_data = EEGout.data;

for ii = 1:nch
    EEGin_data_ch = EEGin_data(ii, :);
    EEGout_data_ch = fftgausshp(EEGin_data_ch, 0.5, fs);
    tmp = EEGout_data_ch(trigs(1) : trigs(end)-1);
    tmp = tmp';
    tmp = reshape(tmp, [], length(trigs)-1);
    tmp = tmp';
    tmp = tmp - avgmat*tmp;
    tmp = tmp';
    tmp = tmp(:)';
    EEGout_data_ch(trigs(1) : trigs(end)-1) = tmp;
    EEGout_data(ii, :) = EEGout_data_ch;
end

EEGout.data = EEGout_data;

end


function avgmat = getavgmat(trigs)

avgmat = zeros(length(trigs)-1);
for ii = 1:length(trigs)-1
    ind = ii-10 : ii+10; ind(ind < 1 | ind > size(avgmat,2)) = [];
    avgmat(ii, ind) = 1/length(ind);
end

end