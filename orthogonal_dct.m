function [U,properties] = orthogonal_dct(n,kind,output)
%ORTHOGONAL_DCT Orthonormal DCT-I--IV and DCT-VIII.
%   U=ORTHOGONAL_DCT(N,KIND) returns a real orthogonal N-by-N matrix.
%   KIND defaults to 2. Type III is the transpose/inverse of type II.
%   Type I requires N >= 2. Type VIII is the odd-period, mixed-boundary
%   symmetric cosine transform with denominator N+1/2, arising also from
%   the Bessel-zero quasi-Hankel kernel at order nu=-1/2. Types V--VII
%   are not implemented. The k-th row is a transform coefficient.
%   Reference: M. Puschel and J. M. F. Moura, IEEE TSP 56 (2008).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [52] https://doi.org/10.1109/TSP.2008.925261
%   Strang, SIAM Review 41 (1999), 135--147.
%   https://doi.org/10.1137/S0036144598336745
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('orthogonal_dct',n,kind);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(kind),kind=2;end
if nargin<3,output='matrix';end
if nargin==2 && (ischar(kind) || isstring(kind)),output=kind;kind=2;end
validateattributes(kind,{'numeric'},{'scalar','integer','>=',1});
if kind~=8
    validateattributes(kind,{'numeric'},{'scalar','integer','<=',4});
end
if kind==1 && n==1
    error('orthogonal_dct:Order','DCT-I requires N >= 2.');
end
if um_is_operator(output)
    U=@(X,flag)um_dct_apply(X,flag,n,kind);
    return
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
    case 8
        U=(2/sqrt(2*n+1))*cos(2*pi*(K+0.5).*(J+0.5)/(2*n+1));
end
end
