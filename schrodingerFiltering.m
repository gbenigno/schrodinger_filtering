function ysf = schrodingerFiltering(y, varargin)
% SCHRODINGERFILTERING Perform Schrödinger filtering on an input vector.
% YSF = SCHRODINGERFILTERING(Y) performs Schrödinger filtering on input
% vector Y and returns output vector YSF.
% 
% YSF = SCHRODINGERFILTERING(Y, FS) performs Schrödinger filtering on input 
% vector Y and returns output vector YSF. FS is the sampling frequency in Hz, 
% and is used to detrend Y by bandpass filtering.
% 
% YSF = SCHRODINGERFILTERING(Y, FS, TRIGS) performs Schrödinger filtering on input
% vector Y and returns output vector YSF. FS is the sampling frequency in Hz, 
% and is used to detrend Y by bandpass filtering. TRIGS 
% is a vector of the indices of the slice acquisition times of Y. The last element of TRIGS 
% is an extra index corresponding to the end of the last slice acquisition, such that the 
% acquisition has a length of LENGTH( TRIGS(1) : (TRIGS(END)-1) ). It is the 
% user's responsibility to ensure that TRIGS adheres to this format. TRIGS is used 
% to partition the slice epochs and apply the function per epoch, which improves its 
% efficiency.
% 
% YSF = SCHRODINGERFILTERING(Y, [], TRIGS) applies Schrödinger filtering on input
% vector Y and returns output vector YSF. The second argument FS--the sampling 
% rate--is left empty, resulting in the detrending step being skipped. TRIGS is a vector
% of the indices of the slice acquisition times of Y. The last element of TRIGS is an extra
% index corresponding to the end of the last slice acquisition, such that the 
% acquisition has a length of LENGTH( TRIGS(1) : (TRIGS(END)-1) ). It is the 
% user's responsibility to ensure that TRIGS adheres to this format. TRIGS is used 
% to partition the slice epochs and apply the function per epoch, which improves its 
% efficiency.
% 
% YSF = SCHRODINGERFILTERING( ___, MF ) also notch-filters the input at the
% mains frequency, MF (Hz), and the first two harmonics.
% 
% NOTE: This function begins to significantly slow for epochs of length ~2000.

% Reference: 

% MIT License
% 
% Copyright (c) 2020 Gabriel Benigno
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% init
[y, lb, ub, t] = sf_init(y);        

% detrend
ydt = sf_detrend(y, varargin);  

% chunk
[ydt_chunk, t_acq] = sf_chunk(ydt, varargin);  
ysf_chunk = sf_chunk(y, varargin);
y_chunk = sf_chunk(y, varargin);

% loop over chunks
for kk = 1 : size(ydt_chunk, 2)
%     fprintf('%u\n',kk)
    
    % rectify
    if iscell(ydt_chunk)
        yp = sf_rec(ydt_chunk{kk}, 1);
        ym = sf_rec(ydt_chunk{kk}, -1);
    else
        yp = sf_rec(ydt_chunk(:,kk), 1);
        ym = sf_rec(ydt_chunk(:,kk), -1);
    end
    t = (1:length(yp))';
    
    % main SF
    if iscell(ysf_chunk)
        ysf_chunk{kk} = sf(yp, ysf_chunk{kk}, 1, lb, ub);
        ysf_chunk{kk} = sf(ym, ysf_chunk{kk}, -1, lb, ub);
    else
        ysf_chunk(:,kk) = sf(yp, ysf_chunk(:,kk), 1, lb, ub);
        ysf_chunk(:,kk) = sf(ym, ysf_chunk(:,kk), -1, lb, ub);
    end
end

% assign output
ysf = y;
if iscell(ysf_chunk)
    ysf(t_acq) = cell2mat(ysf_chunk');
else
    ysf(t_acq) = ysf_chunk(:);
end

end



%-------------------------------------------------------
function ysf = sf(yrec, ysf, pm, lb, ub)

t = 1:length(yrec);

% smooth
[k, yc] = sfsmooth(yrec, lb, ub);

% kernel density estimation
[f, ki] = sfkde(k);

% subtract spikes
while any(islocalmin(f))
    % identify spike components
    [~,locs]=findpeaks(-f, ki);
    spike_components = find(k>locs(1));
    
    % Apodize
    for ii = spike_components(1) : spike_components(end)
        yc(yc(:,ii) < 0.05*max(yc(:,ii)), ii) = 0;
    end   
    
    % Group spike-containing components according to which spikes they represent. From these
    % groups, the domains of the spikes in the input are known. Per spike,
    % null everywhere outside the spike domain so that the spike is
    % reconstructed with maximal accuracy. Then subtract this
    % reconstruction from the input.
    wndws = sf_islands(yc, spike_components);
    
    ctr=0;
    tf=false;
    for ii = 1:size(wndws,2)
        wndw = wndws(:,ii);
        if isequal(wndw, zeros(size(yrec))) || length(wndw(wndw==1)) > 0.5*length(yrec)
            tf=true;
            break
        end
        for ll = 1:4
            try
                h2 = fminbnd(   @(h) mean((SCSA(yrec.*wndw, h) - yrec.*wndw).^2),  lb, ub  );
                MSE = mean((SCSA(yrec.*wndw, h2) - yrec.*wndw).^2);
            catch
                if ll==1
                    ctr=ctr+1;
                end
                continue
            end
            yscsa = SCSA(yrec .* wndw, h2);
            ysf = ysf - pm*yscsa;
            yrec = yrec - yscsa; yrec(yrec < 0)=0;
            if MSE < 1
                break
            end
        end
    end
    if tf
        break
    end
    if ctr == size(wndws,2)
        break
    end
    [k, yc] = sfsmooth(yrec, lb, ub);
    [f, ki] = sfkde(k);
end

end

function [k, yc] = sfsmooth(yrec, lb, ub)
t=(1:length(yrec))';
n=30;
kurts = nan(n,1);
hvec = linspace(lb, ub, n);
for ii = 1:n
    kurts(ii) = kurtosis( SCSA(yrec,hvec(ii)) );
end
h = fminsearch(   @(h) -kurtosis( SCSA(yrec,h) ),  hvec(kurts==max(kurts))  );
[~, ~, ~, k, yc] = SCSA(yrec, h);
end

function [f, ki] = sfkde(k)
[f, ki] =  ksdensity(k);
[~,I] = max(f);
ind = ki<ki(I);
f(ind) = [];
ki(ind) = [];
ind = f>0.8*max(f);
ki(ind)=[];
f(ind)=[];
end

function wndws = sf_islands(yc, spike_components)

islands = sum(yc(:,spike_components),2);
islands(islands>0)=1;
starts = [];
ends=[];
aa=2;
if islands(1) == 1
    starts = [starts; 1];
    for aa = 2 : length(islands)-1
        if islands(aa) == 0
            ends = [ends; aa-1];
            break
        end
    end
end
for ii = aa : length(islands)-1
    if islands(ii) == 1 && islands(ii-1) == 0
        starts = [starts; ii];
        for jj = ii : length(islands)-1
            if islands(jj+1) == 0
                ends = [ends; jj];
                break
            end
        end
    end
end
if islands(end) == 1
    ends = [ends; length(islands)];
end
wndws = zeros(size(yc,1), length(starts));
for ii = 1:length(starts)
    wndws( starts(ii) : ends(ii) ) = 1;
end

end

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

function [y, lb, ub, t] = sf_init(y)
    y=y(:);
    if ~isa(y, 'double')
        y = double(y);
    end
    lb = 0.1;
    ub = 30;    
    t = (1:length(y))';
end

function y_dt = sf_detrend(y, cellarr)

switch length(cellarr)
    
    case 0 % no detrending
        
        y_dt = y;
        
    case {1,2} % detrending by bandpass filtering
        
        fs = cellarr{1};
        if isempty(fs)
            y_dt = y;
        elseif ~isempty(fs)
            y_dt = bandpass(y, [30, 700], fs);
        end
        
    case 3 % detrending by bandpass filtering and notch filtering
        
        fs = cellarr{1};
        mf = cellarr{3};
        if isempty(fs)
            y_dt = y;
        elseif ~isempty(fs)
            y_dt = bandpass(y, [30, 700], fs);
        end
        if ~isempty(mf)
            d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',mf-1,'HalfPowerFrequency2',mf+1, ...
                    'DesignMethod','butter','SampleRate',fs);
            y_dt = filtfilt(d,y_dt);
            d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',2*mf-1,'HalfPowerFrequency2',2*mf+1, ...
                    'DesignMethod','butter','SampleRate',fs);
            y_dt = filtfilt(d,y_dt);
            d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',3*mf-1,'HalfPowerFrequency2',3*mf+1, ...
                    'DesignMethod','butter','SampleRate',fs);
            y_dt = filtfilt(d,y_dt);
        end
        
end

end

function [out, t_acq] = sf_chunk(in, cellarr)

switch length(cellarr)
    case {0,1}
        out = in;
        t_acq = (1 : length(out))';
    case {2,3}
        trigs = cellarr{2};
        if isempty(trigs)
            out = in;
            t_acq = (1 : length(out))';
        elseif ~isempty(trigs)
            L = length(unique(diff(trigs)));
            nep = length(trigs)-1; % number of epochs
            t_acq = trigs(1) : (trigs(end)-1); % indices of acquisition
            if L==1
                out = reshape(in(t_acq), [], nep);
            elseif L > 1
                out = mat2cell(in(t_acq), repelem(diff(trigs),1), 1)';
            end
        end
end

end

function yrec = sf_rec(yin, pm)

yrec = pm*yin;
yrec(yrec<0)=0;

end
