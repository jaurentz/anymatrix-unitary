function [U,properties]=discrete_fourier(n,output)
%DISCRETE_FOURIER Normalized unitary discrete Fourier transform.
%   F(j+1,k+1)=exp(-2*pi*1i*j*k/N)/sqrt(N), j,k=0:N-1.
%   F is symmetric, F^2 reverses indices modulo N, and F^4=I;
%   consequently the exact eigenvalues belong to {1,-1,1i,-1i}.
%   DISCRETE_FOURIER(N,'operator') uses shared FFT multiplication.
%   Reference: P. J. Davis, Circulant Matrices, Wiley, 1979.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [1] https://openlibrary.org/books/OL4407364M/Circulant_matrices
%   [52] https://doi.org/10.1109/TSP.2008.925261
properties={'square','unitary','normal','symmetric','complex', ...
    'scalable','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_dft_apply(X,flag,n);
    return
end
k=(0:n-1)';
U=exp(-2i*pi*(k*k')/n)/sqrt(n);
end
