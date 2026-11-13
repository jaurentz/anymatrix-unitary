function [U,properties] = affine_fourier(n,c1,c2,output)
%AFFINE_FOURIER Unitary chirp-DFT-chirp transform used in AFDM.
%   U=diag(exp(-2i*pi*c2*k^2))*F*diag(exp(-2i*pi*c1*k^2)),
%   k=0:N-1, F(j,k)=exp(-2i*pi*j*k/N)/sqrt(N).
%   C1,C2 are finite real chirp parameters (default zero).
%   This is the modulation transform, not the physical channel matrix.
%   Reference: A. Bemani, N. Ksairi and M. Kountouris, AFDM (2023).
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(c1),c1=0;end
if nargin<3 || isempty(c2),c2=0;end
if nargin<4,output='matrix';end
validateattributes(c1,{'numeric'},{'scalar','real','finite'});
validateattributes(c2,{'numeric'},{'scalar','real','finite'});
k=(0:n-1)';
d1=exp(-2i*pi*c1*k.^2);
d2=exp(-2i*pi*c2*k.^2);
if um_is_operator(output)
    U=@(X,flag)um_affine_apply(X,flag,d1,d2);
    return
end
F=exp(-2i*pi*(k*k')/n)/sqrt(n);
U=bsxfun(@times,d2,bsxfun(@times,F,d1.'));
end
