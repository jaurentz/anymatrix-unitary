function [U,properties]=dpss_transform(n,W)
%DPSS_TRANSFORM Slepian/DPSS orthogonal analysis basis.
%   U=DPSS_TRANSFORM(N,W) returns the rows of the DPSS eigenbasis,
%   ordered by decreasing spectral concentration, using the commuting
%   symmetric tridiagonal of Slepian (1978). 0<W<1/2 is the half-bandwidth.
%   Default W=1/4. Eigenvector signs use maximum-magnitude-positive rule.
%   U is orthogonal. The concentration eigenvalues are NOT the
%   eigenvalues of U, which generally have no simple formula.
%   Only an explicit matrix is returned: a function handle that caches
%   the complete dense eigenbasis would provide no matrix-free benefit.
%   Reference: D. Slepian, Bell Syst. Tech. J. 57 (1978), 1371-1430.
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(W),W=0.25;end
validateattributes(W,{'numeric'},{'scalar','real','finite','>',0,'<',0.5});
j=(0:n-1)';
d=((n-1)/2-j).^2*cos(2*pi*W);
e=(1:n-1)';
e=e.*(n-e)/2;
if n==1
    V=1;
else
    T=diag(d)+diag(e,1)+diag(e,-1);
    [V,D]=eig(T);
    [~,idx]=sort(diag(D),'descend');
    V=V(:,idx);
end
for j=1:n
    [~,pivot]=max(abs(V(:,j)));
    if V(pivot,j)<0,V(:,j)=-V(:,j);end
end
U=V';
end
