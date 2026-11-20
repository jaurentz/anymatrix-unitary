function [U,properties]=krawtchouk_transform(n,p,output)
%KRAWTCHOUK_TRANSFORM Orthonormal finite Krawtchouk polynomial transform.
%   U=KRAWTCHOUK_TRANSFORM(N,P) returns the N-by-N real orthogonal
%   eigenvector basis of the Krawtchouk Jacobi matrix J for degree N-1.
%   0<P<1, default 1/2. Let M=N-1 and j=0,...,M:
%      J(j,j)=M*P+j*(1-2*P),
%      J(j,j+1)=sqrt(P*(1-P)*(j+1)*(M-j)).
%   The eigenvalues of J are the exact grid points 0,...,M.
%   Columns of U are ordered by increasing grid point; the sign of each
%   column is chosen to make its first nonzero component positive.
%   In exact arithmetic the first row is sqrt(w(x)), where
%      w(x)=nchoosek(M,x)*P^x*(1-P)^(M-x).
%   Thus U(j,x) is the orthonormal Krawtchouk function of degree j
%   evaluated at x, with a positive x-dependent sqrt(w(x)) factor.
%
%   IMPORTANT: J has known grid eigenvalues; U itself is orthogonal and
%   has a DIFFERENT eigenvalue distribution. This generator provides an
%   explicit transformation matrix, not J and not its node spectrum.
%   J is diagonalized numerically as a symmetric real matrix. No fast
%   matrix-free Krawtchouk transform or fixed finite trace law is claimed.
%
%   DLMF 18.19 and 18.22: https://dlmf.nist.gov/18.19
%   https://dlmf.nist.gov/18.22
%   Krawtchouk transforms and convolutions, DOI:
%   https://doi.org/10.1007/s13373-018-0132-2
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(p),p=.5;end
validateattributes(p,{'numeric'}, ...
    {'scalar','real','finite'});
if ~(p>0 && p<1)
    error('krawtchouk_transform:Probability','P must lie strictly between 0 and 1.');
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('krawtchouk_transform:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('krawtchouk_transform',n);
    U.reference='Krawtchouk transforms and convolutions, doi:10.1007/s13373-018-0132-2';
    U.notes=['The associated symmetric Jacobi matrix has exact ' ...
        'eigenvalues 0:N-1 with binomial spectral weights at e1. ' ...
        'These are not the eigenvalues or trace weights of the ' ...
        'orthogonal transformation matrix returned by this generator.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('krawtchouk_transform:Mode', ...
        'Only explicit matrix and spectral modes are implemented.');
end
if n==1,U=1;return;end
M=n-1;
j=(0:M)';
diagonal=M*p+j*(1-2*p);
a=sqrt(p*(1-p)*(1:M)'.*(M-(1:M)'+1));
J=diag(diagonal)+diag(a,1)+diag(a,-1);
[U,D]=eig(J);
[~,order]=sort(real(diag(D)));
U=real(U(:,order));
for k=1:n
    if U(1,k)<0,U(:,k)=-U(:,k);end
end
end
