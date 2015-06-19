%{
    Test for deconvolution of Wigner estimates.
%}

addpath('toolbox/', 'SlantStackP/');

%% 
% Parameters.

% sampling size
N = 256;
% the observation domain is [-r,r] x [r,r]
r = 5;
% noise level
eta = .85;
eta = .9;
% gamma in [0,1/4]
gamma = (1-eta)/(4*eta);

%%
% Helpers.

radon = @(f)real( FastSlantStack(f) );
% iradon = @(R)real( Inv_FastSlantStack(R) );
mynorm = @(x)norm(x(:));
normalize = @(R)R ./ repmat( sum(R), [size(R,1) 1] );

name = 'vacuum';
name = 'schrodinger-cat';
name = 'single-photon';

rep0 = ['results/' name '/'];
if not(exist(rep0))
    mkdir(rep0);
end

%%
% Load target Wigner function.

rho = load_rho(name, N, r);
% normalize through radon transforms
R = radon(rho);
R = normalize(max(R,0));
rho = iradon(R);

% display
for mt = {'surf' 'ls'}
    clf;
    rho_surf(rho, mt{1});  drawnow; 
    saveas(gcf, [rep0 'rho-' mt{1} '.png'], 'png');
end
clf;

%% 
% Total #samples. Controls amount of noise. 

% C. BUTUCEA paper:
% one-photon : n=5000 eta=.9
% shrod, n=500000 eta=.95 / eta=.85

n = 100e6; % 
n = 5000; % 1-photon
n = 500000;  % shrod
n = 100000;  % 1-photon
str = ['eta' num2str(round(eta*100)) '-' 'n' num2str(round(n/1e3)) 'e3'];

rep = [rep0 '/' str '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Blurring level (measure in sampling locations).
% [-2*r,2*r] for [-N,N] samples noise level in real unit is sigma * 2*r/N

sigma = sqrt(2*gamma) * N / (2*r); 

%%
% Do the sampling.

sampling_mode = 'samples';
sampling_mode = 'poisson';
R1 = perform_radon_sampling(R, n, sigma, sampling_mode);
rho1 = iradon(R1);

if 0
clf;
rho_surf(rho1); drawnow;
saveas(gcf, [rep 'rho-' str '-noisy-surf.png'], 'png');
end

%%
% Load deconvolution/denoising operator

window_type = 'hard';
window_type = 'smooth';
iRad = load_deconvolution_operator(N, @iradon, window_type);

% kappa should be rather compared with a gaussian of width 2*N/sigma or so
% where sigma=sqrt( 2*gamma ) * N / (2*r)

param_mode = 'kappa';
param_mode = 'h';


M0 = ceil( sqrt(log(n/(2*gamma))) ); % theoretical number
M = 20; % number of tested values

kappa_list = [];
switch name
    case 'schrodinger-cat'
        switch str
            case 'eta95-n500e3'
                kappa_list = linspace(20,40,M); %  n=.5e6, eta=.95
            case 'eta90-n500e3'
                kappa_list = linspace(16,35,M); %  n=.5e6, eta=.9 HARD
                kappa_list = linspace(12,22,M); %  n=.5e6, eta=.9 SOFT
            case 'eta85-n500e3'
                kappa_list = linspace(17,30,M); %  n=.5e6, eta=.85
                kappa_list = linspace(12,20,M); %  n=.5e6, eta=.85 SOFT
        end
    case 'single-photon'
        switch str
            case 'eta90-n1e3'
                kappa_list = linspace(5,20,M); %  n=5000, eta=.9
            case 'eta90-n100e3'
                kappa_list = linspace(8,17,M); %  n=100000, eta=.9
        end
end
hlist = (N/(2*r)) ./ kappa_list;

if strcmp(param_mode, 'h')
    hlist = linspace(min(hlist),max(hlist), M);
    kappa_list = (N/(2*r)) ./ hlist;
end

if 0
    M = 10;
    hlist = linspace(3.2, 1.5, M); % 1/2 * (1 -(0:M-1)'/M);
    kappa_list = (N/(2*r)) ./ hlist;
end

myerr = @(rho1)mynorm(rho - rho1)/mynorm(rho);
func = @(kappa)myerr(iRad(R1,sigma,kappa) );

Rho = [];
for i=1:M
    progressbar(i,M);
    rho1 = iRad(R1,sigma,kappa_list(i));
    Rho(:,i) = rho1(:); 
end

% error in Linf and l^2
ErrLinf = max( abs(Rho-repmat(rho(:), [1 M])) ) / max(abs(rho(:)));
ErrL2   = sqrt( sum( (Rho-repmat(rho(:), [1 M])).^2 ) ) / mynorm(rho);

% this is the kappa of the article.
k = 1;
x = log(M0); %   c~1
rn = @(x)max( sqrt((1+x)/n), (1+x)/n );

% Lepski matrix
normi = @(x)2*N*norm(x(:), 'inf'); % normalization needed
LkMat = zeros(M);
LkMat0 = zeros(M);
LkMat1 = zeros(M);
for m=1:M
    for j=m:M
        LkMat0(m,j) = normi(Rho(:,m) - Rho(:,j));
        LkMat1(m,j) = - 2*k*rn(x+log(M)) * ( exp(gamma/hlist(j)^2)-exp(gamma/hlist(m)^2) );
    end
end
LkMat = LkMat0 + k*LkMat1;
% Lepski functional
Lk = max(LkMat, [], 2);
% optimal choice
[~,m] = min(Lk);

myplot = @(x,y,c)plot(x,y,c,'LineWidth', 2);
y = ErrLinf;
switch param_mode
    case 'kappa'
        x = kappa_list;
    case 'h'
        x = hlist;
end

figure(2);
clf; hold on;
myplot(x, y, 'k');
axis tight; 
% axis([min(x) max(x) min(y) 1]);
box on;
set(gca, 'FontSize', 15);
saveas(gcf, [rep 'error-trial-' str '.eps'], 'epsc');

% oracle result
[~,i] = min(ErrLinf); kappa = kappa_list(i);
rhoA = iRad(R1,sigma,kappa);
figure(1);
for mt = {'surf' 'ls'}
    clf;
    rho_surf(rhoA, mt{1}); drawnow;
    saveas(gcf, [rep 'rho-' str '-oracle-' mt{1}  '.png'], 'png');
end
clf;
