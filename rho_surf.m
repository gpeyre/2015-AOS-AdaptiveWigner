function rho_surf(rho, drawmode, options)

% rho_surf - display Wigner function in 3D
%
%   rho_surf(rho, drawmode, options);
%
%   drawmode is either 'surf' or 'levelset'
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<2
    drawmode = 'surf';
end

options.null = 0;
q = getoptions(options, 'nbls', 12); % #levelset displayed

% normalize in [-1,1]
rho = rho/max(abs(rho(:)));
% rho(rho==min(rho(:)))=-1;
% rho(rho==max(rho(:)))=+1;
        
switch drawmode
    case 'surf'
        h = surf(rho); % , 'CDataMapping','direct'
        shading interp;
        colormap parula(256);
        caxis([-1 1]);
        axis off;
        camlight;
    case {'levelset' 'ls'}
        N = size(rho,1);
        t = linspace(-1,1,N+1); t(end)=[];
        u = linspace(-1,1,q+2); u([1 end]) = [];
        hold on;
        imagesc(t,t, rho);
        colormap parula(256);        
        caxis([-1 1]);
        contour(t,t,rho, u, 'k', 'LineWidth', 2);
        axis image; axis off;
    otherwise
        error('Unknown mode');
end

end