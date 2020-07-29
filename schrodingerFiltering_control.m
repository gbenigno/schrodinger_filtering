function ysf = schrodingerFiltering_control(y)
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
y=y(:);
t=(1:length(y))';
lb=0.1; ub=30;
ysf=y;

% rectify
yp = y; yp(yp<0) = 0; yp_orig = yp;
ym = -y; ym(ym<0)=0; ym_orig = ym;
% spks_p = spks; spks_p(spks_p<0) = 0;
% spks_m = -spks; spks_m(spks_m<0) = 0;

% smooth (positive)
h = fminbnd(   @(h) -kurtosis( SCSA(yp,h) ),  lb, ub  );
[~, ~, ~, k, ycp] = SCSA(yp, h);

% kernel density estimation (positive)
[f, ki] = ksdensity(k);

% subtract spikes (positive)
while any(islocalmin(f)) && any(islocalmin(ksdensity(kurtosis(ycp))))
    [~,locs] = findpeaks(-f, ki);
    spike_components = find(k>locs(1));
    
    % find if any secondary components in spike_components
    note=[];
    for ii = 1:length(spike_components)
        [~,locs] = findpeaks(ycp(:,ii),'minpeakheight',0.01*max(ycp(:,ii)));
        if length(locs) == 2 
            [~,I] = max(ycp(:,spike_components));
            pair = [find(locs(1) < I & locs(2) > I); ii];
            note = [note, pair];
            wndw = zeros(size(yp));
            tmp = sum(ycp(:,pair),2);
            wndw(tmp>0.01*max(tmp)) = 1;
            h2 = fminbnd(   @(h) mean((SCSA(yp.*wndw,h) - yp.*wndw).^2),  lb, ub  );
            yscsa = SCSA(yp.*wndw, h2);
            ysf = ysf - yscsa;
            yp = yp - yscsa; yp(yp<0)=0;
        end
    end
    
    ind = setdiff(spike_components, note(:));
    for ii = ind(1) : ind(end)
        wndw = zeros(size(yp));
        wndw(ycp(:,ii)>0.01*max(ycp(:,ii))) = 1;
        h2 = fminbnd(   @(h) mean((SCSA(yp.*wndw,h) - yp.*wndw).^2),  lb, ub  );
        yscsa = SCSA(yp.*wndw, h2);
        ysf = ysf - yscsa;
        yp = yp - yscsa; yp(yp<0)=0;
    end

    h = fminbnd(   @(h) -kurtosis( SCSA(yp,h) ),  lb, ub  );
    [~, ~, ~, k, ycp] = SCSA(yp, h);
    [f, ki] = ksdensity(k);
end

% smooth (negative)
h = fminbnd(   @(h) -kurtosis( SCSA(ym,h) ),  lb, ub  );
[~, ~, ~, k, ycm] = SCSA(ym, h);

% kernel density estimation (negative)
[f, ki] = ksdensity(k);

% subtract spikes (negative)
while any(islocalmin(f)) && any(islocalmin(ksdensity(kurtosis(ycm))))
    [~,locs] = findpeaks(-f, ki);
    spike_components = find(k>locs(1));
    
    % find if any secondary components in spike_components
    note=[];
    for ii = 1:length(spike_components)
        [~,locs] = findpeaks(ycm(:,ii),'minpeakheight',0.01*max(ycm(:,ii)));
        if length(locs) == 2 
            [~,I] = max(ycm(:,spike_components));
            pair = [find(locs(1) < I & locs(2) > I); ii];
            note = [note, pair];
            wndw = zeros(size(ym));
            tmp = sum(ycm(:,pair),2);
            wndw(tmp>0.01*max(tmp)) = 1;
            h2 = fminbnd(   @(h) mean((SCSA(ym.*wndw,h) - ym.*wndw).^2),  lb, ub  );
            yscsa = SCSA(ym.*wndw, h2);
            ysf = ysf + yscsa;
            ym = ym - yscsa; ym(ym<0)=0;
        end
    end
    
    ind = setdiff(spike_components, note(:));
    for ii = ind(1) : ind(end)
        wndw = zeros(size(ym));
        wndw(ycm(:,ii)>0.01*max(ycm(:,ii))) = 1;
        h2 = fminbnd(   @(h) mean((SCSA(ym.*wndw,h) - ym.*wndw).^2),  lb, ub  );
        yscsa = SCSA(ym.*wndw, h2);
        ysf = ysf + yscsa;
        ym = ym - yscsa; ym(ym<0)=0;
    end

    h = fminbnd(   @(h) -kurtosis( SCSA(ym,h) ),  lb, ub  );
    [~, ~, ~, k, ycm] = SCSA(ym, h);
    [f, ki] = ksdensity(k);
end

end





%-------------------------------------------------------
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
try
    I = simp(psi_n.^2, dx); % Normalization of the eigenfunction 
catch
    save('~/Muller Lab Dropbox/gabriel/schrodinger_filtering/simulation/spike_simulation/wkspc.mat')
end
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
