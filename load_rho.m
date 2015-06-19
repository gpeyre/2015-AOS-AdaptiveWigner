function [rho,R] = load_rho(name, N, r)

% load_rho - load a Wigner function
%
% rho = load_rho(name, N, R0)
%
%   rho is the Wigner function
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<3
    r = 4; % radius of observation
end

t = linspace(-r,r,N+1); t(end)=[];
[Y,X] = meshgrid(t,t);
gauss = @(x,y)exp(-(X-x).^2-(Y-y).^2);
A = 2*X.^2+2*Y.^2;
switch name
    case 'vacuum'
        rho = gauss(0,0);
    case 'single-photon'
        rho = -gauss(0,0) .* (1-A);
    case 'two-photon'
        rho = gauss(0,0) .* (A.^2 + 4*A + 2); % or +/-4*A ...
    case 'schrodinger-cat'
        q0 = 3;
        rho = 1/2 * gauss(q0,0) + 1/2 * gauss(-q0,0) + ...
            gauss(0,0) .* cos(2*q0*Y);
    otherwise
        error('Unknown');
end

end