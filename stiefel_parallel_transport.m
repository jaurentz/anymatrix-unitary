function [Q,properties]=stiefel_parallel_transport(n,p,A,B,t,mode)
%STIEFEL_PARALLEL_TRANSPORT Intrinsic canonical Stiefel tangent isometry.
%   Q=STIEFEL_PARALLEL_TRANSPORT(N,P,A,B,T) returns the real
%   D-by-D orthogonal matrix of Levi-Civita parallel transport along
%   a geodesic on St(p,n), D=N*P-P*(P+1)/2. This is NOT the ambient
%   N-by-N frame lift STIEFEL_GEODESIC.
%
%   The base frame is Y0=[I_p;0], with horizontal so(n) generator
%       K=[A,-B'; B,zeros(n-p)],
%   where A is real skew-symmetric p-by-p and B is (n-p)-by-p.
%   The geodesic is Y(t)=expm(t*K)*Y0 (canonical quotient metric).
%   At Y(t), take tangent basis vectors expm(t*K)*F_ij*Y0,
%   where F_ij=e_i*e_j'-e_j*e_i', 1<=i<=p, i<j<=n.
%   These are orthonormal for g_Y(V,W)=tr(V'*(I-Y*Y'/2)*W).
%
%   For the normal homogeneous quotient SO(n)/SO(n-p), the connection
%   in the moving horizontal frame is the projected Lie bracket:
%      C_(ij,ab)=([K,F_ab])_(i,j),  S=-C/2,  Q=expm(t*S).
%   C is real skew-symmetric; hence Q is EXACTLY orthogonal in
%   mathematical arithmetic. This is intrinsic tangent transport.
%   The geodesic velocity is parallel. n=p specializes to SO(n);
%   p=1 is spherical geodesic parallel transport (Q=I in the
%   specified moving frame).
%
%   MODE='matrix' (default) returns Q; MODE='model' returns Q, K, S,
%   orthonormal horizontal basis index pairs, Y0 and Yt. MODE='spectral'
%   labels general finite transform eigenphases UNKNOWN (no closed law).
%   Dense matrix/model: explicit D-by-D connection and expm, O(D^2)
%   structure constants, O(D^3) exponential and O(D^2) storage.
%   Require N<=64,D<=256 for these dense reference modes.
%
%   MODE='operator': matrix-free scaled degree-24 Taylor action of
%   exp(t*S), using ONLY compact (A,B), without forming either n-by-n
%   K or d-by-d S. Every S times vector costs O(N*P^2) and O(N*P)
%   work storage (fixed P<=32). For ||S||<=||K||_F, splitting into
%   m=max(1,ceil(abs(t)*||K||_F)) steps gives a Taylor remainder
%   at most m*exp(1)/25!*(1+exp(1)/25!)^(m-1) in real arithmetic.
%   This APPROXIMATE action is not exactly orthogonal; unlike
%   the dense mathematical expm it incurs polynomial and roundoff
%   errors. Forward/adjoint/transpose/multi-RHS all supported.
%   Operator allows N<=1024, D<=16384, P<=32, m<=4096.
%   This is an explicit projected-Lie-bracket method with a
%   Taylor action; it does NOT claim to reproduce Nguyen--Sommer's
%   more sophisticated large-time algorithm or improve its bound.
%
%   [97] Edelman, Arias, Smith (1998),
%   https://doi.org/10.1137/S0895479895290954
%   Nguyen, Sommer (2025), https://doi.org/10.1137/24M1684189
properties={'square','unitary','normal','orthogonal','real', ...
    'parameter-dependent'};
if nargin==0,Q=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
    &&n==floor(n)&&n>=2&&n<=1024)
    error('stiefel_parallel_transport:Order', ...
        'N must be an integer from 2 to 1024.');
end
if nargin<2||isempty(p),p=min(2,n);end
if ~(isnumeric(p)&&isscalar(p)&&isreal(p)&&isfinite(p) ...
    &&p==floor(p)&&p>=1&&p<=n)
    error('stiefel_parallel_transport:Order', ...
        'P must be an integer with 1<=P<=N.');
end
d=n*p-p*(p+1)/2;
if nargin<3||isempty(A)
    A=zeros(p);
    if p>1,A(1,2)=-.27;A(2,1)=.27;end
end
if ~isnumeric(A)||~isreal(A)||~isequal(size(A),[p,p]) ...
    ||any(~isfinite(A(:))) ...
    ||norm(A+A','fro')>1e-11*max(1,norm(A,'fro'))
    error('stiefel_parallel_transport:Skew', ...
        'A must be finite real skew-symmetric P-by-P.');
end
if nargin<4||isempty(B)
    B=zeros(n-p,p);
    if n>p
        B(1,1)=.32;
        if p>1,B(1,2)=.17;end
    end
end
if ~isnumeric(B)||~isreal(B)||~isequal(size(B),[n-p,p]) ...
    ||any(~isfinite(B(:)))
    error('stiefel_parallel_transport:Tangent', ...
        'B must be finite real (N-P)-by-P.');
end
if nargin<5||isempty(t),t=1;end
if ~(isnumeric(t)&&isscalar(t)&&isreal(t)&&isfinite(t))
    error('stiefel_parallel_transport:Time', ...
        'T must be a finite real scalar.');
end
if nargin<6||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('stiefel_parallel_transport:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','spectral','operator'}))
    error('stiefel_parallel_transport:Mode', ...
        'MODE must be matrix, model, spectral, or operator.');
end
if any(strcmp(mode,{'matrix','model'})) && (n>64||d>256)
    error('stiefel_parallel_transport:Order', ...
        'Dense reference requires N<=64 and D<=256.');
end
if strcmp(mode,'operator') && (p>32||d>16384)
    error('stiefel_parallel_transport:OperatorSize', ...
        'Operator requires P<=32 and D<=16384.');
end
if strcmp(mode,'spectral')
    Q=unitary_spectral_reference('stiefel_parallel_transport',d);
    Q.reference='Edelman-Arias-Smith (1998); Nguyen-Sommer (2025)';
    Q.notes=['Canonical Stiefel tangent parallel transport, D-by-D; ', ...
       'generic finite trace eigenphase law unknown. ', ...
       'Do not confuse geodesic-generator eigenvalues with Q eigenphases.'];
    return
end

ii=zeros(d,1);jj=zeros(d,1);
h=0;
for i=1:p
    for j=i+1:n
        h=h+1;ii(h)=i;jj(h)=j;
    end
end
if strcmp(mode,'operator')
    Q=@(X,flag)um_stiefel_connection_expmv( ...
        X,flag,n,p,A,B,t,ii,jj);
    return
end
K=[A,-B';B,zeros(n-p)];
% For F_ab=e_a e_b'-e_b e_a', take projected bracket coordinates
% [K,F_ab](i,j). In this canonical horizontal orthonormal basis,
% -ad_K/2 is the Levi-Civita parallel transport generator.
S=zeros(d,d);
for col=1:d
    a=ii(col);b=jj(col);
    z=K(ii,a).*(jj==b)-K(ii,b).*(jj==a) ...
      -(ii==a).*reshape(K(b,jj),[],1) ...
      +(ii==b).*reshape(K(a,jj),[],1);
    S(:,col)=-.5*z;
end
S=(S-S')/2; % enforce exact skew matrix within floating roundoff
U=expm(t*S);
if strcmp(mode,'matrix')
    Q=U;
else
    G=expm(t*K);
    Q=struct('U',U,'connection_generator',S, ...
        'horizontal_geodesic_generator',K, ...
        'basis_row',ii,'basis_col',jj,'dimension',d, ...
        'Y0',eye(n,p),'Yt',G(:,1:p), ...
        'moving_frame',G,'time',t,'metric','canonical', ...
        'notes',['U acts on orthonormal moving tangent coordinates. ', ...
         'It is intrinsic parallel transport, not an ambient lift. ', ...
         'Dense prototype O(D^3) setup; fast exponential action ', ...
         'and noncanonical metrics are not implemented.']);
end
end
