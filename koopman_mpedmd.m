function [Q,properties]=koopman_mpedmd(n,PX,PY,weights,mode)
%KOOPMAN_MPEDMD Gram-whitened measure-preserving EDMD unitary.
%   Q=KOOPMAN_MPEDMD(N,PX,PY,WEIGHTS,MODE) is a STRICTLY unitary
%   N-by-N data-derived Koopman approximation from the Colbrook (2023)
%   mpEDMD orthogonal-Procrustes construction. PX and PY are M-by-N
%   dictionary EVALUATIONS at paired snapshots x_j and y_j=T(x_j),
%   respectively: PX(j,k)=psi_k(x_j), PY(j,k)=psi_k(y_j).
%   The user may build PX/PY from arbitrary sample states/dictionaries.
%   WEIGHTS is a positive M-vector (default uniform).
%
%   The default data use M=max(64,2*N+7) equispaced circle snapshots
%   x_j=2*pi*j/M, rotation y_j=x_j+0.37, and complex Fourier
%   dictionary psi_k(x)=exp(1i*(k-1)*x), k=1..N.
%   In the default case the finite Fourier dictionary is invariant
%   under this circle rotation; its exact finite phases are
%       exp(1i*(k-1)*0.37), k=1..N.
%   This is a reproducible FINITE truncation, not a generic proof of
%   infinite Koopman spectrum from arbitrary data.
%
%   With W=diag(weights/sum(weights)), set
%       G=PX'*W*PX (Hermitian SPD)
%       A=PX'*W*PY (cross Gram)
%       B=G^(-1/2)*A*G^(-1/2)
%       Q=polar(B)=U*V'  where B=U*S*V' (full SVD).
%   This Q is unitary in Euclidean coordinates, exactly the database
%   matrix. The associated COEFFICIENT-SPACE Koopman model
%       K=G^(-1/2)*Q*G^(1/2)
%   satisfies K'*G*K=G but is generally NOT Euclidean unitary.
%   Ordinary EDMD K_edmd=G\A is not automatically G-unitary.
%
%   This orientation is verified against the author's MATLAB
%   mpEDMD(G,A): it takes SVD(G^(-1/2)*A'*G^(-1/2)) and builds
%   V*U', equivalent to polar(B) in the above convention.
%
%   G must be strictly SPD and cond(G)<=1e11, validated.
%   B may be rank-deficient: full SVD U*V' supplies a unitary
%   completion, but its action on null directions is NONUNIQUE.
%   The 'metric' output exposes numerical rank and completion flag.
%
%   MODE='matrix' (default): return the Euclidean unitary Q.
%   MODE='metric': return G,A,B,Q,K,EDMD, H=G^(1/2), Hinv,
%       W, singular values, cross rank and nonuniqueness flag.
%   MODE='data': return snapshot dictionaries PX/PY, normalized
%       weights, default snapshot angles (when applicable) and
%       the corresponding metric model.
%   MODE='edmd': return G\A, NOT necessarily unitary.
%   MODE='spectral': EXACT finite eigenphases/moments for the
%       built-in invariant Fourier/rotation example ONLY; else
%       UNKNOWN analytical finite phase law, despite the ability
%       to compute numerical eig(Q). No universal Koopman density.
%
%   An invertible measure-preserving transformation induces a
%   unitary infinite-dimensional Koopman operator. A noninvertible
%   measure-preserving transformation instead induces an isometry
%   not necessarily onto, and an arbitrary finite EDMD projection
%   need not be a unitary. Colbrook's mpEDMD enforces a finite
%   metric isometry by the data-derived Procrustes problem.
%
%   [148] M. J. Colbrook, SIAM J. Numer. Anal. 61 (2023), 1585-1608:
%       https://doi.org/10.1137/22M1521407
%   Reference source implementation:
%       https://github.com/MColbrook/Measure-preserving-Extended-Dynamic-Mode-Decomposition
%   [146] Colbrook and Townsend, CPAM 77 (2024), 221-283:
%       https://doi.org/10.1002/cpa.22125
properties={'square','unitary','normal','complex', ...
            'parameter-dependent'};
if nargin==0,Q=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
        &&n==floor(n)&&n>=1&&n<=64)
    error('koopman_mpedmd:Order','N must be 1..64.');
end
defaultData=(nargin<2||isempty(PX))&&(nargin<3||isempty(PY)) ...
     &&(nargin<4||isempty(weights));
if nargin<2,PX=[];end
if nargin<3,PY=[];end
if nargin<4,weights=[];end
if nargin<5||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('koopman_mpedmd:Mode','MODE must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','metric','data','spectral','edmd'}))
    error('koopman_mpedmd:Mode','Unknown output mode.');
end
if isempty(PX)&&isempty(PY)
    count=max(64,2*n+7);
    x=2*pi*(0:count-1)'/count;
    theta=.37;
    y=mod(x+theta,2*pi);
    k=0:n-1;
    PX=exp(1i*x*k);
    PY=exp(1i*y*k);
else
    if isempty(PX)||isempty(PY)
        error('koopman_mpedmd:Snapshots','PX and PY must be supplied together.');
    end
    x=[];y=[];theta=[];
end
if ~(isnumeric(PX)&&isnumeric(PY)&&isfloat(PX)&&isfloat(PY) ...
    &&ismatrix(PX)&&ismatrix(PY)&&isequal(size(PX),size(PY)) ...
    &&size(PX,2)==n&&size(PX,1)>=n ...
    &&size(PX,1)<=20000 ...
    &&all(isfinite(PX(:)))&&all(isfinite(PY(:))))
    error('koopman_mpedmd:Dictionary', ...
        'PX, PY must be finite M-by-N dictionary evaluations, M>=N.');
end
count=size(PX,1);
if isempty(weights),weights=ones(count,1)/count;end
if ~(isnumeric(weights)&&isvector(weights)&&isreal(weights) ...
    &&numel(weights)==count&&all(isfinite(weights(:))) ...
    &&all(weights(:)>0))
    error('koopman_mpedmd:Weights', ...
        'Require one finite strictly positive weight per snapshot.');
end
weights=double(weights(:));
weights=weights/sum(weights);
% Compute measured inner products using actual positive quadrature.
WPX=bsxfun(@times,weights,PX);
WPY=bsxfun(@times,weights,PY);
G=PX'*WPX;
G=(G+G')/2;
A=PX'*WPY;
[Vg,Dg]=eig(G,'vector');
ev=real(Dg);
if min(ev)<=1e-11*max(ev) || any(~isfinite(ev)) ...
        || max(ev)/min(ev)>1e11
    error('koopman_mpedmd:Gram', ...
        'Dictionary weighted Gram must be well-conditioned SPD.');
end
H=(Vg.*sqrt(ev).')*Vg';
Hinv=(Vg.*(1./sqrt(ev)).')*Vg';
B=Hinv*A*Hinv;
[Ub,S,Vb]=svd(B);
Q=Ub*Vb';
rankCross=sum(diag(S)>1e-11*max(1,max(diag(S))));
if strcmp(mode,'matrix'),return;end
K=Hinv*Q*H;
EDMD=G\A;
if strcmp(mode,'edmd')
    Q=EDMD;
    return
end
metric=struct('G',G,'A',A,'B',B,'Q',Q,'K',K, ...
    'EDMD',EDMD,'H',H,'Hinv',Hinv, ...
    'weights',weights,'singular_values',diag(S), ...
    'cross_rank',rankCross,'rank_deficient_completion',rankCross<n, ...
    'normalization','normalized positive sample weights', ...
    'notes',['Q is Euclidean unitary; K is G-unitary, ', ...
       'whereas EDMD need not be. For singular B, full SVD U*V'' ', ...
       'is a nonunique null-space unitary completion.']);
if strcmp(mode,'metric')
    Q=metric;return
end
if strcmp(mode,'data')
    Q=struct('PX',PX,'PY',PY,'weights',weights, ...
        'snapshot_x',x,'snapshot_y',y,'theta',theta, ...
        'dictionary','user-specified evaluations or default Fourier characters', ...
        'metric',metric);
    return
end
blank=struct('kind','unknown','density',[],'support',[], ...
    'atoms',[],'moments',[]);
finite=struct('kind','unknown','eigenvalues',[], ...
    'moments',[],'atoms',[]);
notes=['Arbitrary data/finite dictionary gives an EXACTLY UNITARY ', ...
    'Procrustes matrix Q but no asserted analytic finite eigenphase ', ...
    'measure of an infinite Koopman operator. EDMD itself need ', ...
    'not be G-unitary; metric K is.'];
if defaultData
    z=exp(1i*.37*(0:n-1)');
    finite.kind='exact';
    finite.eigenvalues=@()z;
    finite.moments=@(k)reshape(mean(bsxfun(@power,z, ...
       reshape(k,1,[])),1),size(k));
    notes=['For the default exact equispaced rotation with invariant ', ...
       'Fourier dictionary, the finite Koopman phases are ', ...
       'exp(+i*k*.37), k=0..N-1, in the pullback convention. ', ...
       'This does not establish general infinite Koopman spectral ', ...
       'convergence from noisy data.'];
end
Q=struct('family','koopman_mpedmd','dimension',n, ...
    'normalization','probability (density with respect to dtheta)', ...
    'finite',finite,'limit',blank,'ensemble',blank,'underlying',blank, ...
    'reference','Colbrook, mpEDMD, SIAM J. Numer. Anal. 61 (2023)', ...
    'notes',notes,'cross_rank',rankCross);
end
