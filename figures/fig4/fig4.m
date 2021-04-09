clearvars;clc

addpath ../../external/linspecer

%% spiky signal
TRsl = 0.1; pausepercent = 0.2; stf = 700; nch = 1; Nsl = 0; fs=2048;
[t, ~, ~, ~, Ns] = gradsim(fs, TRsl, Nsl, pausepercent, stf); % gradient artifact
wno = dsp.ColoredNoise('white',Ns,nch);
x = wno();
spks = zeros(size(x));
bfr = round(0.05*length(t));
jj=1;
while jj < 5
    loc = randi([bfr, length(t)-bfr+1]);
    spktmp = 30*rms(x)*gaussmf(t, [1/(2048), t(loc)]);
    if ~isempty(intersect(find(spks>0.01*30*rms(x)), find(spktmp>0.01*30*rms(x))))
        continue
    end
    spks = spks + spktmp;
    jj=jj+1;
end
jj=1;
while jj < 5
    loc = randi([bfr, length(t)-bfr+1]);
    spktmp = 30*rms(x)*gaussmf(t, [1/(2048), t(loc)]);
    if ~isempty(intersect(find(abs(spks)>0.01*30*rms(x)), find(abs(spktmp)>0.01*30*rms(x))))
        continue
    end
    spks = spks - spktmp;
    jj=jj+1;
end
y = x+spks;

figure
plot(t,y,'k','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% flipped
figure
plot(t,-y,'k','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% rectify positive
figure
yp = y;
yp(yp<0)=0;
plot(t,yp,'k','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% rectify negative
figure
ym = -y;
ym(ym<0)=0;
plot(t,ym,'k','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% smooth yp
lb=0.1; ub=30;
n=30;
kurts = nan(n,1);
hvec = linspace(lb, ub, n);
for ii = 1:length(hvec)
    kurts(ii) = kurtosis( SCSA(yp,hvec(ii)) );
end
h = fminsearch(   @(h) -kurtosis( SCSA(yp,h) ),  hvec(kurts==max(kurts))  );
[~, ~, ~, k, ycp] = SCSA(yp, h);

figure
plot(t,yp,'k',t,SCSA(yp,h),'r','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% smoothed yp spectrum prob dens fn
figure
[f,ki]=ksdensity(k);
plot(ki,f,'k','linewidth',2)
set(gca,'TickLength',[0 0])
set(gca,'fontsize',16)
xlabel('eigenvalue')
ylabel('probability density')
xlim([ki(1) ki(end)])
[~,loc]=findpeaks(-f,ki,'NPeaks',1);
hold on
area(ki(ki>loc),f(ki>loc),'FaceColor','r','FaceAlpha',0.3)
xline(loc,'--')
hold off

%% smoothed yp spectrum
figure
C=linspecer(length(find(k>loc)));
for ii = 1:length(find(k>loc))
    stem(ii,k(ii),'filled','color',C(ii,:),'linewidth',2)
    hold on
end
for ii = length(find(k>loc))+1 : length(k)
    stem(ii,k(ii),'filled','color','k','linewidth',2)
end
xlim([0 length(k)+1])
yline(loc,'--')
xlabel('component number')
ylabel('eigenvalue')
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'TickLength',[0 0])
box on
hold off

%% selected components of smoothed yp
figure
kp = k;
locp=loc;
plot(t,yp,'k',t,ycp(:,kp>loc),'linewidth',2)
colororder(C)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])
set(gca,'ytick',[],'xtick',[])

%% smooth ym
figure
lb=0.1; ub=30;
n=30;
kurts = nan(n,1);
hvec = linspace(lb, ub, n);
for ii = 1:length(hvec)
    kurts(ii) = kurtosis( SCSA(ym,hvec(ii)) );
end
h = fminsearch(   @(h) -kurtosis( SCSA(ym,h) ),  hvec(kurts==max(kurts))  );
[~, ~, ~, k, ycm] = SCSA(ym, h);
plot(t,ym,'k',t,SCSA(ym,h),'r','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])

%% smoothed ym spectrum prob dens fn
figure
[f,ki]=ksdensity(k);
plot(ki,f,'k','linewidth',2)
set(gca,'TickLength',[0 0])
set(gca,'fontsize',16)
xlabel('eigenvalue')
ylabel('probability density')
xlim([ki(1) ki(end)])
[~,loc]=findpeaks(-f,ki,'NPeaks',1);
hold on
area(ki(ki>loc),f(ki>loc),'FaceColor','r','FaceAlpha',0.3)
xline(loc,'--')
hold off

%% smoothed ym spectrum
figure
C=linspecer(length(find(k>loc)));
for ii = 1:length(find(k>loc))
    stem(ii,k(ii),'filled','color',C(ii,:),'linewidth',2)
    hold on
end
for ii = length(find(k>loc))+1 : length(k)
    stem(ii,k(ii),'filled','color','k','linewidth',2)
end
xlim([0 length(k)+1])
yline(loc,'--')
xlabel('component number')
ylabel('eigenvalue')
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'TickLength',[0 0])
box on
hold off

%% selected components of smoothed yp
figure
km = k;
locm=loc;
plot(t,ym,'k',t,ycm(:,km>loc),'linewidth',2)
colororder(C)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])
set(gca,'ytick',[],'xtick',[])

%% flip ym comps
figure
plot(t,-ym,'k',t,-ycm(:,km>loc),'linewidth',2)
colororder(C)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])
set(gca,'ytick',[],'xtick',[])

%% subtract from y
figure
cd ../..
ysf = schrodingerFiltering(y);
plot(t,y,'k',t,ysf,'r','linewidth',2)
set(gca,'ytick',[],'xtick',[])
xlim([t(1) t(end)])
set(gca,'ytick',[],'xtick',[])


%% --------------------------------------------------

function [y_scsa, Nh, psi_n_nor, kappa, y_comps] = SCSA(y, h, dx, gm)
%%% courtesy of Dr. Meriem Taous-Kiriati
% y is the signal of interest
% dx is the spacing of the x-values of y(x) (should be 1)
% gm is lowercase gamma from eq 8 of article (should be 0.5)

if nargin==2
    dx=1;
    gm=0.5;
end

Lcl = gamma(gm+1) / ( 2*sqrt(pi)*gamma(gm+3/2) ); % classical constant (eq 8 in article)

% remove the negative part
ymin = min(y); 
y_nn = y - ymin; % non-negative

% Build second-order differentiation matrix 
D2=diffmat2o(y_nn, dx);

% construct and solve the eigenvalue problem
H = -h*h*D2 - diag(y_nn); % The Schrodinger operator
[psi, lambda] = eig(H); % solve
all_lambda = diag(lambda); % generate vector of eigenvalues from the given diagonal matrix of eigenvalues
ind = find(all_lambda<0); 
neg_lambda = all_lambda(ind); % find which entries are negative eigenvalues
kappa = diag( (abs(neg_lambda)).^gm ) ; % get kappa values by finding the sqrt of the abs value of the negative eigenvalues
Nh = length(kappa); % number of negative eigenvalues
psi_n = psi(:,ind(:,1)); % The eigenfunctions of the negative eigenvalues
I = simp(psi_n.^2, dx); % Normalization of the eigenfunction 
psi_n_nor = psi_n./sqrt(I);  % The L^2 normalized eigenfunction 
y_comps = (h/Lcl)*(  (psi_n_nor.^2) * kappa    ) .^ ( 2 / (1+2*gm) ) ;
y_scsa = sum( y_comps, 2 ) + ymin;

if size(y) ~= size(y_scsa)
    y_scsa = y_scsa';
end

kappa = diag(kappa);

end








function D2=diffmat2o(y, Delta_x)
%**********************************************************************
%*********             second-order differentiation matrix           *********
%**********************************************************************
% 
%Author: Zineb Kaisserli
% modified by Gabriel Benigno to be consistent with notation in SCSA paper (Math. Control Signals Syst. (2013) 25:37?61);
% also got rid of Fs argument; made =1 always
M = length(y);
Delta = 2*pi/M;
ex = kron( (M-1):-1:1, ones(M,1) ); % 'j-k' from eqns 35 and 36 of article
if mod(M,2)==0
    dx = -pi^2/(3*Delta^2)-(1/6)*ones(M,1);
    test_bx = -(-1).^ex*(0.5)./(sin(ex*Delta*0.5).^2);
    test_tx =  -(-1).^(-ex)*(0.5)./(sin((-ex)*Delta*0.5).^2);
else
    dx = -pi^2/(3*Delta^2)-(1/12)*ones(M,1);
    test_bx = -0.5*((-1).^ex).*cot(ex*Delta*0.5)./(sin(ex*Delta*0.5));
    test_tx = -0.5*((-1).^(-ex)).*cot((-ex)*Delta*0.5)./(sin((-ex)*Delta*0.5));
end
Ex = full(spdiags([test_bx dx test_tx],[-(M-1):0 (M-1):-1:1],M,M));
D2=(Delta/Delta_x)^2*Ex;

end








function y = simp(f,dx)
%**********************************************************************
%*********              Numerical integration                 *********
%**********************************************************************
% 
% Author: Taous Meriem Laleg
%  This function returns the numerical integration of a function f
%  using the Simpson method

n=length(f);
I(1)=1/3*f(1)*dx;
I(2)=1/3*(f(1)+f(2))*dx;

for i=3:n
    if(mod(i,2)==0)
        I(i)=I(i-1)+(1/3*f(i)+1/3*f(i-1))*dx;
    else
        I(i)=I(i-1)+(1/3*f(i)+f(i-1))*dx;
    end
end
y=I(n);
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

rmpath ../../external/linspecer
