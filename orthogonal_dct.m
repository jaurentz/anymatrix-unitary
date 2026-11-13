function [U,properties] = orthogonal_dct(n,kind)
%ORTHOGONAL_DCT Orthonormal discrete cosine transform of type I, II, III, IV.
%   U=ORTHOGONAL_DCT(N,KIND) returns a real orthogonal N-by-N matrix.
%   KIND defaults to 2. Type III is the transpose/inverse of type II.
%   Type I requires N >= 2. The k-th row is a transform coefficient.
%   Reference: M. Puschel and J. M. F. Moura, IEEE TSP 56 (2008).
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(kind),kind=2;end
validateattributes(kind,{'numeric'},{'scalar','integer','>=',1,'<=',4});
if kind==1 && n==1
    error('orthogonal_dct:Order','DCT-I requires N >= 2.');
end
[K,J]=ndgrid(0:n-1,0:n-1);
switch kind
    case 1
        c=ones(n,1);c([1 n])=1/sqrt(2);
        U=sqrt(2/(n-1))*(c*c').*cos(pi*K.*J/(n-1));
    case {2,3}
        c=ones(n,1);c(1)=1/sqrt(2);
        T=sqrt(2/n)*bsxfun(@times,c,cos(pi*K.*(J+0.5)/n));
        if kind==2,U=T;else,U=T';end
    case 4
        U=sqrt(2/n)*cos(pi*(K+0.5).*(J+0.5)/n);
end
end
