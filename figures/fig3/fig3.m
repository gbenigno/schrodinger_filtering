clearvars;clc
%%
% SCSA decomp of gaussian
figure
t = 15:0.025:35;
y = 10*gaussmf(t,[0.5 25]);
h0=1;
fun = @(h) mean( (SCSA(y, h) - y).^2 );
h = fminsearch(fun, h0);
[~, Nh, ~, kappa, yc] = SCSA(y, h);
C=linspecer(Nh);
plot(t,y,'k','linewidth',2)
hold on
plot(t,yc)
colororder(C)
set(gca,'XTick',[],'ytick',[])
hold off

%% closeup of decomp
figure
C=linspecer(Nh);
plot(t,y,'k')
hold on
plot(t,yc)
colororder(C)
xlim([23,27])
ylim([0,2.7])
set(gca,'XTick',[],'ytick',[])
hold off

%% first 4 comps
figure
C=linspecer(4);
plot(t,y,'k')
hold on
plot(t,yc(:,1:4),'linewidth',2)
colororder(C)
width=1000;
height=700;
xlim([24.5,25.5])
ylim([0,2.7])
set(gca,'XTick',[],'ytick',[])
hold off

%% eigenspectrum
figure
C=linspecer(Nh);
for ii=1:Nh
    stem(ii,kappa(ii),'filled','color',C(ii,:),'linewidth',2)
    hold on
end
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'ytick',[],'TickLength',[0 0])
hold off

%% noisy input and smooth reconstruction
figure
t = 1:0.05:50;
gc=10*gaussmf(t,[0.5 25]);
yc = randn(size(t)) + gc;

h0=1;
fun = @(h) mean( (SCSA(yc, h) - gc).^2 );
h = fminsearch(fun, h0);

[ys, ~, ~, kappa, ~] = SCSA(yc, h);

p=plot(nan,nan,'k',nan,nan,'r',t,yc,'k',t,ys,'r');
p(1).LineWidth=8;
p(2).LineWidth=8;
p(3).LineWidth=2;
p(4).LineWidth=2;
set(gca,'fontsize',14)
set(gca,'ytick',[],'xtick',[])
legend('noisy input',['SCSA reconstruction (h=',num2str(round(h)),')'])

%% eigenspectrum
figure
stem(kappa,'filled','color','r','linewidth',2)
xlim([0 length(kappa)+1])
xlabel('component number')
ylabel('eigenvalue')
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'ytick',[],'TickLength',[0 0])


%% noisy input and full reconstruction
figure
h0=1;
fun = @(h) mean( (SCSA(yc, h) - yc).^2 );
h = fminsearch(fun, h0);

[ys, ~, ~, kappa, ~] = SCSA(yc, h);

p=plot(nan,nan,'k',nan,nan,'r',t,yc,'k',t,ys,'r');
p(1).LineWidth=8;
p(2).LineWidth=8;
p(3).LineWidth=2;
p(4).LineWidth=2;
set(gca,'fontsize',14)
set(gca,'ytick',[],'xtick',[])
legend('noisy input',['SCSA reconstruction (h=',num2str(round(h)),')'])

%% eigenspectrum
figure
stem(kappa,'filled','color','r','linewidth',2)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'ytick',[],'TickLength',[0 0])

%% very noisy input (3 signal peaks)
figure
t = 1:0.05:70;
m=2;
w=0.3;
g1=m*gaussmf(t,[w 25]);
g2=m*gaussmf(t,[w 10]);
g3=m*gaussmf(t,[w 50]);
gb=g1+g2+g3;

yb = randn(size(t)) + gb;

p=plot(t,yb,'k');
p(1).LineWidth = 1;
pbaspect([3 1 1])
set(gca,'ytick',[],'xtick',[])

%% noisy input (3 peaks) and first 3 schr comps and true signal peaks
figure
h0=1;
fun = @(h) mean( (SCSA(yb, h) - gb).^2 );
h = fminsearch(fun, h0);

[ys, Nh, psi_n_nor, kappa, ycomps] = SCSA(yb, h/7);

p=plot(t,yb,'k',t,ycomps(:,1:3),'r',t,gb,'g');
p(1).LineWidth = 1;
p(2).LineWidth = 2;
p(3).LineWidth=2;
p(4).LineWidth=2;
p(5).LineWidth=2;

pbaspect([3 1 1])
set(gca,'ytick',[],'xtick',[])
set(gca,'fontsize',14)

%% eigenspectrum
figure
stem(kappa,'filled','color','r','linewidth',1)
xlim([0, length(kappa)+1])
xlabel('component number')
ylabel('eigenvalue')
set(gca,'fontsize',16)
xlabel('component number')
ylabel('eigenvalue')
set(gca,'ytick',[],'TickLength',[0 0])
pbaspect([3 1 1])



%--------------------------
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