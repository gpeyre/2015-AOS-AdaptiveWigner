Toolbox to reproduce the numerical results of the paper

Minimax and adaptive estimation of the Wigner function in quantum homodyne tomography with noisy data
K. Lounici, K. Meziani and G. Peyré
Preprint 2015

* test_wigner.m : main script that performs sampling, deblurring/denoising and parameter selection.
* compute_rho.m : 
* load_deconvolution_operator.m : load the linear estimator
* load_rho.m : load a Wigner function
perform_radon_sampling.m
* rho_surf.m : display of a Wigner function.
* SlantStackP/ : Fast Slant Stack Transform toolbox implementing
    Fast Slant Stack: A notion of Radon Transform for Data in a Cartesian Grid which is Rapidly Computible, Algebraically Exact, Geometrically Faithful and Invertible
    A. Averbuch, R.R. Coifman, D.L. Donoho, M. Israeli, J. Walden
which makes use of the polar FFT detailed in
    Fast and accurate Polar Fourier transform
    A. Averbuch, R.R. Coifman, D.L. Donoho, M. Elad, M. Israeli
    Applied and Computational Harmonic Analysis, Volume 21, Issue 2, September 2006, Pages 145-167
