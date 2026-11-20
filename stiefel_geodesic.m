function [G,properties]=stiefel_geodesic(n,p,A,R,frames,t,output)
%STIEFEL_GEODESIC Ambient unitary lift of canonical Stiefel geodesic.
%   G=STIEFEL_GEODESIC(N,P,A,R,FRAMES,T) requires 2P<=N and
%   FRAMES=[Y,Q] with 2P orthonormal columns. A is P-by-P
%   skew-Hermitian, and R is any P-by-P matrix. The tangent
%   H=Y*A+Q*R generates canonical-metric Stiefel geodesic
%   Y(t)=[Y,Q]*exp(t*[A,-R';R,zeros(P)])*[eye(P);zeros(P)].
%   The N-by-N ambient unitary is identity on span(FRAMES)^perp.
%   Spectrum: exp(T*eig(K)) plus N-2P ones, with the eigenproblem
%   solved only for small K of order 2P. These spectral values are
%   mathematically exact but numerical in the generic case.
%   The ambient lift is not intrinsic Stiefel tangent transport.
%   Operator mode applies the precomputed 2P block in O(NP).
%   [93] Edelman, Arias, Smith, SIAM J. Matrix Anal. Appl. (1998).
%   https://doi.org/10.1137/S0895479895290954
properties={'square','unitary','normal','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,G=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(p),p=min(2,floor(n/2));end
validateattributes(p,{'numeric'},{'scalar','integer','positive','finite'});
if 2*p>n,error('stiefel_geodesic:Order','N must be at least 2P.');end
if nargin<3||isempty(A)
    A=zeros(p);
    if p>1,A(1,2)=-.3;A(2,1)=.3;end
end
if ~isnumeric(A)||~isequal(size(A),[p,p])||any(~isfinite(A(:))) ...
        ||norm(A+A','fro')>1e-11*max(1,norm(A,'fro'))
    error('stiefel_geodesic:Skew','A must be a skew-Hermitian P-by-P matrix.');
end
if nargin<4||isempty(R),R=diag(.3+.2*(1:p));end
if ~isnumeric(R)||~isequal(size(R),[p,p])||any(~isfinite(R(:)))
    error('stiefel_geodesic:Tangent','R must be a finite P-by-P matrix.');
end
if nargin<5||isempty(frames),frames=eye(n,2*p);end
if ~isnumeric(frames)||~isequal(size(frames),[n,2*p]) ...
        ||any(~isfinite(frames(:))) ...
        ||norm(frames'*frames-eye(2*p),'fro')>1e-10*max(1,p)
    error('stiefel_geodesic:Frames','FRAMES must be N-by-2P orthonormal.');
end
if nargin<6||isempty(t),t=1;end
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<7||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('stiefel_geodesic:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
K=[A,-R';R,zeros(p)];
if strcmp(output,'spectral')
    z=[exp(t*eig(K));ones(n-2*p,1)];
    G=unitary_spectral_reference('stiefel_geodesic',n);
    G.finite.kind='exact';G.finite.eigenvalues=@()z;
    G.finite.moments=@(k)reshape(arrayfun(@(j) ...
        mean(z.^j),k(:)),size(k));
    G.reference='Edelman-Arias-Smith (1998), doi:10.1137/S0895479895290954';
    G.notes='Exact spectral reduction to 2P-by-2P skew-Hermitian generator; small eig used numerically.';
elseif any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    B=expm(t*K);
    if any(strcmp(output,{'operator','implicit'}))
        G=@(X,flag)um_stiefel_apply(X,flag,frames,B);
    else
        G=eye(n)+frames*(B-eye(2*p))*frames';
    end
else
    error('stiefel_geodesic:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
