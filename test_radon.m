%{
    test for Radon transform
%}

mynorm = @(x)norm(x(:));

name = 'parrot-1';
name = 'square';

N = 256;
addpath('toolbox/', 'SlantStackP/');
addpath('radon/');
% Inv_FastSlantStack / FastSlantStack

f = load_image(name, N);
f = rescale(sum(f,3));

% f = zeros(N); f(end/2-2:end/2+2,:) = 1;

% FSS: -pi/4 -> 3*pi/4 equaly spaced in tangent
tic;
warning off;
R = real( FastSlantStack(f) );
f1 = real( Inv_FastSlantStack(R) );
warning on;
fprintf('FSS: Time=%.3f\n', toc);
fprintf('FSS: Error=%.3f\n', mynorm(f-f1)/mynorm(f));


Theta = linspace(0,180,2*N+1); Theta(end) = [];
tic;
R1 = radon(f,Theta);
f2 = iradon(R1,Theta, 'Spline', 'cosine'); %, 'nearest', 'Hann' );
f2 = clamp(f2(2:end-1,2:end-1));
fprintf('Rad: Time=%.3f\n', toc);
fprintf('Rad: Error=%.3f\n', mynorm(f-f2)/mynorm(f));
