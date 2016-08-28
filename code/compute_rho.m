function rho = compute_rho(name, n, t)

% compute_rho - compute rho matrix
%
%   rho = compute_rho(name, n, t);
%   
%   t is the (optional) parameter of the distribution
%
%   Copyright (c) 2012 Gabriel Peyre

[K,J] = meshgrid(0:n-1,0:n-1);
f = @(i)factorial(i);

switch name
    case 'vacuum'
        rho = zeros(n); rho(1,1) = 1;
    case 'single-photon'
        rho = zeros(n); rho(2,2) = 1;
    case 'coherent'
        q0 = t/sqrt(2);
        rho = exp(-abs(q0)^2) * q0.^(J+K) ./ sqrt( f(J).*f(K) );
    case 'squeezed'
        error('Not yet implemented');
    case 'thermal'
        beta = t;
        rho = (1-exp(-beta)) * exp( -beta*K ) .* (K==J);
    case 'shrodinger-cat'
        q0 = t;%  check formula bellow
        a = 2*(q0/sqrt(2)).^(J+K) ./ sqrt(f(K).*f(J)) / (exp(q0^2/2)+exp(-q0^2/2));
        rho = (mod(J,2)==0) .* (mod(K,2)==0) .* a;
    otherwise 
        error('Unknown');
end

end