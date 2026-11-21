function [Q,properties]=real_flag_parallel_transport(blockSizes,K,t,mode)
%REAL_FLAG_PARALLEL_TRANSPORT Intrinsic normal-metric real flag isometry.
%   Q=REAL_FLAG_PARALLEL_TRANSPORT(BLOCKS,K,T) is the D-by-D real
%   orthogonal Levi-Civita parallel transport along a geodesic on
%   F=SO(n)/S(O(k1) x ... x O(kr)), sum(BLOCKS)=n, r>=2.
%   The tangent dimension is D=sum_{a<b}ka*kb and the normal-homogeneous
%   metric is inherited from <X,Y>=tr(X'*Y)/2 on so(n).
%
%   The moving orthonormal tangent frame is indexed by all skew
%   F_ij=e_i e_j'-e_j e_i' with 1<=i<j<=n in DIFFERENT BLOCKS.
%   K is a real skew-Hermitian horizontal n-by-n matrix: each diagonal
%   block of K must be zero. The flag geodesic is exp(t*K)*H and the
%   LC connection is S=-project_m(ad_K)/2 in this moving frame:
%       S_[(i,j),(a,b)] = -[K,F_ab]_(i,j)/2.
%   Then Q=expm(t*S), mathematically Q'*Q=I_D.
%   In the 2-block (real Grassmann) symmetric-space case [m,m] is
%   vertical, S=0 and Q=I in the geodesic moving frame. In the
%   full-flag all-one blocks case, m=so(n) and the expression
%   specializes to Lie-group bi-invariant parallel transport.
%
%   MODE='matrix' (default) gives Q. MODE='model' gives Q, the
%   projected connection generator S, horizontal basis index pairs,
%   the initial/final orthogonal projectors for all flag subspaces,
%   the moving frame, and K. MODE='spectral' gives an exact identity
%   spectrum for 2-block flags and UNKNOWN general eigenphases
%   otherwise. Neither geodesic generator eigenvalues nor principal
%   angles are misidentified as Q eigenphases.
%
%   This is an explicit dense prototype, 3<=n<=32, D<=256.
%   S assembly O(D^2), exponential O(D^3), Q storage/action O(D^2).
%   Not a Nguyen-Sommer low-rank matrix exponential-action algorithm,
%   not the embedded non-normal metric, and not complex flags.
%
%   [168] Nguyen and Sommer, SIAM J. Matrix Anal. Appl. (2025)
%   https://doi.org/10.1137/24M1684189
properties={'square','unitary','orthogonal','real','normal', ...
    'parameter-dependent'};
if nargin==0,Q=[];return;end
if ~(isnumeric(blockSizes)&&isreal(blockSizes) ...
    &&isvector(blockSizes)&&numel(blockSizes)>=2 ...
    &&all(isfinite(blockSizes(:)))&&all(blockSizes(:)>=1) ...
    &&all(blockSizes(:)==fix(blockSizes(:))))
    error('real_flag_parallel_transport:Blocks', ...
        'BLOCKS must be a vector of at least two positive integers.');
end
blockSizes=double(blockSizes(:)');
n=sum(blockSizes);
if n<3||n>32
    error('real_flag_parallel_transport:Order', ...
        'The prototype requires 3<=sum(BLOCKS)<=32.');
end
r=numel(blockSizes);
d=(n^2-sum(blockSizes.^2))/2;
if d<1||d>256
    error('real_flag_parallel_transport:Order', ...
        'Tangent dimension must be at most 256.');
end
labels=zeros(n,1);
pos=0;
for a=1:r
    labels(pos+(1:blockSizes(a)))=a;
    pos=pos+blockSizes(a);
end
pairI=zeros(d,1);pairJ=zeros(d,1);h=0;
for i=1:n
    for j=i+1:n
        if labels(i)~=labels(j)
            h=h+1;pairI(h)=i;pairJ(h)=j;
        end
    end
end
if h~=d,error('real_flag_parallel_transport:Internal', ...
        'Unexpected horizontal dimension.');end
if nargin<2||isempty(K)
    K=zeros(n);
    for j=1:d
        a=pairI(j);b=pairJ(j);
        K(a,b)=.38*sin(.31*j)+.17*cos(.27*j);
        K(b,a)=-K(a,b);
    end
end
if ~(isnumeric(K)&&isreal(K)&&isequal(size(K),[n n]) ...
    &&all(isfinite(K(:))) ...
    &&norm(K+K','fro')<=1e-11*max(1,norm(K,'fro')))
    error('real_flag_parallel_transport:Skew', ...
        'K must be finite, real, skew, n-by-n.');
end
for a=1:r
    ids=find(labels==a);
    if norm(K(ids,ids),'fro')>1e-11*max(1,norm(K,'fro'))
        error('real_flag_parallel_transport:Horizontal', ...
            'K must have zero within-block diagonal skew components.');
    end
end
if nargin<3||isempty(t),t=1;end
if ~(isnumeric(t)&&isscalar(t)&&isreal(t)&&isfinite(t))
    error('real_flag_parallel_transport:Time', ...
        'T must be finite real scalar.');
end
if nargin<4||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('real_flag_parallel_transport:Mode', ...
        'MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','spectral'}))
    error('real_flag_parallel_transport:Mode', ...
        'Only matrix, model, or spectral is supported.');
end
if strcmp(mode,'spectral')
    Q=unitary_spectral_reference('real_flag_parallel_transport',d);
    Q.block_sizes=blockSizes;
    Q.reference='Nguyen--Sommer 2025, DOI 10.1137/24M1684189';
    if r==2
        Q.finite.kind='exact';
        Q.finite.eigenvalues=@()ones(d,1);
        Q.finite.moments=@(k)ones(size(k));
        Q.notes=['Two-block Grassmann symmetric-space case: ', ...
            'normal-metric geodesic moving-frame parallel ', ...
            'transport is exactly identity.'];
    else
        Q.notes=['Normal-homogeneous real flag tangent transport; ', ...
            'generic finite orthogonal eigenphase law unknown. ', ...
            'Geodesic K eigenvalues are unrelated to Q phases.'];
    end
    return
end
S=zeros(d);
for col=1:d
    a=pairI(col);b=pairJ(col);
    z=K(pairI,a).*(pairJ==b)-K(pairI,b).*(pairJ==a) ...
      -(pairI==a).*reshape(K(b,pairJ),[],1) ...
      +(pairI==b).*reshape(K(a,pairJ),[],1);
    S(:,col)=-z/2;
end
S=(S-S')/2;
Q=expm(t*S);
if strcmp(mode,'model')
    G=expm(t*K);
    P0=cell(r,1);Pt=cell(r,1);
    for a=1:r
        ids=find(labels==a);
        P=zeros(n);P(ids,ids)=eye(numel(ids));
        P0{a}=P;Pt{a}=G*P*G';
    end
    Q=struct('U',Q,'connection_generator',S, ...
        'horizontal_geodesic_generator',K, ...
        'basis_row',pairI,'basis_col',pairJ, ...
        'dimension',d,'block_sizes',blockSizes, ...
        'initial_projectors',{P0},'endpoint_projectors',{Pt}, ...
        'moving_frame',G,'time',t,'metric','normal homogeneous', ...
        'notes',['Intrinsic real flag LC tangent parallel transport ', ...
          'in moving orthonormal quotient frame. Dense prototype ', ...
          'O(D^3) exponential; NOT ambient flag rotation, ', ...
          'embedded-metric transport or fast matrix-free action.']);
end
end
