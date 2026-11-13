function [U,properties] = fresnel_transform(n,output)
%FRESNEL_TRANSFORM Normalized circulant discrete Fresnel transform.
%   The first column is exp(-1i*pi/4)*exp(1i*pi*(k+delta)^2/n)/sqrt(n),
%   k=0:N-1, where delta=0 for even N and 1/2 for odd N. Circular
%   translations give the other columns. The parity correction makes the
%   chirp N-periodic. The exact eigenvalues are FFT of the first column.
%   Application: orthogonal chirp-division multiplexing (OCDM).
%   Reference: X. Ouyang and J. Zhao, IEEE Trans. Commun. 64 (2016).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [49] https://doi.org/10.1109/TCOMM.2016.2594792
properties={'square','unitary','normal','circulant','complex', ...
    'scalable','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
k=(0:n-1)';
delta=mod(n,2)/2;
c=exp(-1i*pi/4+1i*pi*(k+delta).^2/n)/sqrt(n);
if nargin<2,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_fft_apply(X,flag,'circulant',fft(c),[]);
    return
end
idx=mod((0:n-1)'-(0:n-1),n)+1;
U=c(idx);
end
