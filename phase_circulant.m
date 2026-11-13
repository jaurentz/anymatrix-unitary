function [U,properties] = phase_circulant(n,theta,output)
%PHASE_CIRCULANT Circulant unitary with prescribed eigenphases.
%   U=PHASE_CIRCULANT(N,THETA) has eigenvalues exp(1i*THETA(j)), where
%   THETA is a real N-vector in Fourier frequency order. In the unitary
%   Fourier convention, U=F'*diag(exp(1i*THETA))*F.
%   Default THETA(j)=2*pi*(j-1)/N, giving a circulant shift.
%   Reference: P. J. Davis, Circulant Matrices, Wiley, 1979.
properties={'square','unitary','normal','circulant','complex', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0, U=[]; return; end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(theta), theta=2*pi*(0:n-1)'/n; end
validateattributes(theta,{'numeric'},{'vector','real','finite'});
if numel(theta)~=n, error('phase_circulant:Size','THETA must have N entries.'); end
if nargin<3,output='matrix';end
lambda=exp(1i*theta(:));
if um_is_operator(output)
    U=@(X,flag)um_fft_apply(X,flag,'circulant',lambda,[]);
    return
end
column=ifft(lambda);
indices=mod((0:n-1)'-(0:n-1),n)+1;
U=column(indices);
end
