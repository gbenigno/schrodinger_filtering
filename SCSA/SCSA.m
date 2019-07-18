%%% courtesy of Dr. Meriem Taous-Kiriati of King Abdullah University of Science and Technology

function [y_scsa, Nh, psi_n_nor, kappa, y_comps] = SCSA(y, h, dx, gm)

% y is the signal of interest
% dx is the spacing of the x-values of y(x) (should be 1)
% gm is lowercase gamma from eq 8 of article (should be 0.5)

%%
Lcl = gamma(gm+1) / ( 2*sqrt(pi)*gamma(gm+3/2) ); % classical constant (eq 8 in article)

%% remove the negative part
ymin = min(y); 
y_nn = y - ymin; % non-negative

%% Build second-order differentiation matrix 
D2=diffmat2o(y_nn, dx);

%% construct and solve the eigenvalue problem
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