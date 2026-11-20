function [Q,properties]=lie_geodesic_transport(n,kappa,t,output)
%LIE_GEODESIC_TRANSPORT Left-trivialized transport along a geodesic on U(n).
%   For K=1i*diag(KAPPA), gamma(t)=expm(t*K) is a geodesic of the
%   bi-invariant metric of U(n). In left-trivialized coordinates parallel
%   transport on its real Lie algebra u(n) is exp(-t*ad_K/2).
%   Basis: n diagonal i*E_jj followed by each off-diagonal pair
%   (E_jk-E_kj)/sqrt(2), 1i*(E_jk+E_kj)/sqrt(2), j<k.
%   The resulting N^2-by-N^2 matrix is real sparse orthogonal, with
%   spectrum 1 (N copies), exp(+/-1i*t*(kappa_j-kappa_k)/2).
%   Normalized Fourier moments: |sum(exp(1i*l*t*kappa/2))|^2/N^2.
%   This is left-trivialized transport, NOT tangent mapping in fixed
%   ambient matrix coordinates.
%   Bi-invariant canonical connection reference [95]:
%   https://doi.org/10.1109/LCSYS.2024.3405485
properties={'square','unitary','normal','orthogonal','real','sparse', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,Q=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(kappa),kappa=.35*(0:n-1)';end
if ~isnumeric(kappa)||~isreal(kappa)||~isvector(kappa)|| ...
        numel(kappa)~=n||any(~isfinite(kappa(:)))
    error('lie_geodesic_transport:Kappa','KAPPA must be N real finite entries.');
end
kappa=kappa(:);
if nargin<3||isempty(t),t=1;end
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('lie_geodesic_transport:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
gaps=zeros(n*(n-1)/2,1);z=0;
for j=1:n-1
    for k=j+1:n
        z=z+1;gaps(z)=kappa(j)-kappa(k);
    end
end
alpha=t*gaps/2;
if strcmp(output,'spectral')
    ev=[ones(n,1);exp(1i*alpha);exp(-1i*alpha)];
    Q=unitary_spectral_reference('lie_geodesic_transport',n);
    Q.dimension=n*n;Q.finite.kind='exact';Q.finite.eigenvalues=@()ev;
    Q.finite.moments=@(ell)reshape(arrayfun(@(k) ...
        abs(sum(exp(.5i*k*t*kappa)))^2/n^2,ell(:)),size(ell));
    Q.reference='Canonical Lie connection; doi:10.1109/LCSYS.2024.3405485';
    Q.notes='Exact eigenphase-difference law; left-trivialized real orthogonal transport.';
elseif any(strcmp(output,{'operator','implicit'}))
    Q=@(X,flag)um_lie_transport_apply(X,flag,n,alpha);
elseif any(strcmp(output,{'matrix','explicit'}))
    M=numel(alpha);a=n+2*(0:M-1)'+1;b=a+1;
    c=cos(alpha);s=sin(alpha);
    rows=[(1:n)';a;a;b;b];
    cols=[(1:n)';a;b;a;b];
    Q=sparse(rows,cols,[ones(n,1);c;s;-s;c],n*n,n*n);
else
    error('lie_geodesic_transport:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
