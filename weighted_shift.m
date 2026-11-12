function [U,properties] = weighted_shift(n,theta)
%WEIGHTED_SHIFT Sparse phase-weighted cyclic shift (monomial unitary).
%   U*e_j=exp(1i*THETA(j))*e_(j+1), with indices modulo N.
%   The eigenvalues are the Nth roots of exp(1i*sum(THETA)), and
%   U^N=exp(1i*sum(THETA))*I. Default THETA(j)=pi*(j-1)/N.
%   Reference: P. J. Davis, Circulant Matrices, Wiley, 1979.
properties={'square','unitary','normal','sparse','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0, U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(theta), theta=pi*(0:n-1)'/n;end
validateattributes(theta,{'numeric'},{'vector','real','finite'});
if numel(theta)~=n, error('weighted_shift:Size','THETA must have N entries.');end
j=1:n;
U=sparse(mod(j,n)+1,j,exp(1i*theta(:)).',n,n);
end
