function R1 = perform_radon_sampling(R, n, sigma, sampling_mode)

% perform_radon_sampling - computes sampled Radon histograms
%
%   R1 = perform_radon_sampling(R, n, sigma, sampling_mode);
%
%   R is the noiseless map.
%   R1 is the noisy map computed by drawing n/size(R,2) samples of each 
%       R(:,i), that are supposed to be normalized density (i.e.
%       should sum to 1) and adding a randn*sigma additive noise.
%   Thus R1 is a convolved and noisy version of R.
%   sampling_mode is either 'samples' where the sampled is really performed
%   and 'poisson' where the sampling is simulated by computing a Poisson
%   noise. 
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<4
    sampling_mode = 'samples';
end

N = size(R,1)/2;
Q = round(n/(2*N));  % #samples/orientation

normalize =@(R)R ./ repmat( sum(R), [size(R,1) 1] );
% fourier filtering
filt = @(R,h)real(ifft( fft(R) .* fft(h) ) );
t = [0:N-1 -N:-1]';
replicate = @(h)h(:)*ones(1,2*N);
gauss1d = @(sigma)normalize(exp(-t.^2/(2*sigma^2)));        
gauss2d = @(sigma)replicate( gauss1d(sigma) );       

switch sampling_mode
    case 'samples'
        %% real sampling %%
        % create using sampling        
        I0 = randdisc(R,Q);
        % add Gaussian noise and then quantize
        I1 = I0 + randn(size(I0))*sigma;
        I1 = clamp(round(I1),1, 2*N);
        % histograms
        R0 = R * 0; R1 = R0;
        for k=1:size(I0,2)
            R0(:,k) = hist( I0(:,k), 1:2*N )'/Q;
            R1(:,k) = hist( I1(:,k), 1:2*N )'/Q;
        end
    case 'poisson'
        %% poisson approximation -- faster sampling %%
        gfilt = @(R,sigma)filt(R,gauss2d(sigma));
        % Add the convolution.
        R0 = gfilt(R, sigma);
        % generate Poisson noise
        lambda = max(1e-9,R0*Q);
        R1 = poissrnd(lambda)/Q;  
    otherwise
        error('Unknown sampling mode');
end

end