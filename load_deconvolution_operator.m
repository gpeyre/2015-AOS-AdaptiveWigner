function [iRad,ifilt] = load_deconvolution_operator(N, iradon, window_type)

% load_deconvolution_operator - load deconvolution operators
%
%   [iRad,ifilt] = load_deconvolution_operator(N, iradon, window_type);
%
%   N is the width of the image
%   iRad(R,sigma,kappa) compute the deconvolved-regularized radon transform,
%       where sigma is the blurring width (assumed to be a Gaussian
%       convolution) and kappa is the denoising parameter (the larger, the
%       more denoising), which corresponds to the maximal keep 1D Fourier
%       frequency. 
%
%   Copyright (c) 2015 Gabriel Peyre

t = [0:N-1 -N:-1]';
switch window_type
    case 'hard'
        w = @(kappa)double(abs(t)<=kappa); % windowing
    case 'smooth'
        maxnorm = @(h)h/max(h(:));
        wf = @(kappa,t) (t<=kappa) + (t>kappa & t<=2*kappa) .* cos( pi/2*(t-kappa)/kappa ).^2;
        w = @(kappa)wf(kappa,abs(t));
    otherwise
        error('Unknown window type');
end
filt = @(R,h)real(ifft( fft(R) .* fft(h) ) );
replicate = @(h)h(:)*ones(1,2*N);
normalize =@(R)R ./ repmat( sum(R), [size(R,1) 1] );
myifft = @(h)( real( ifft(h) ) );
myinv = @(h)1 ./ max(real(h),1e-10); 
gauss1d = @(sigma)normalize(exp(-t.^2/(2*sigma^2)));  
igauss = @(sigma,kappa)myifft( w(kappa) .*  myinv(fft(gauss1d(sigma))) );
ifilt = @(R,sigma,kappa)filt(R, replicate(igauss(sigma,kappa)) );
iRad = @(R,sigma,kappa)iradon(ifilt(R, sigma,kappa ));

end