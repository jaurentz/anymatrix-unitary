function [U,properties] = block_cmv_pencil(n,q,alpha,beta,mode)
%BLOCK_CMV_PENCIL Strict unitary block-CMV pencil, order N*Q.
%   U=BLOCK_CMV_PENCIL(N,Q,ALPHA,BETA) uses an array ALPHA of shape
%   Q-by-Q-by-(N-1) of strictly contractive matrix Verblunsky blocks
%   and a terminal unitary BETA (Q-by-Q). Default ALPHA=0, BETA=I.
%   L,M are block-diagonal unitaries with C=L*M, and
%   P = FEVAL(NAME,...,'pencil') returns A=M and B=L' so that B\A=C.
%   The 2Q-by-2Q local unitary uses
%       [alpha',sqrt(I-alpha'*alpha);
%        sqrt(I-alpha*alpha'),-alpha].
%   The terminal block is beta'. The eigenvalues lie on the circle.
%   For Q=1 the result agrees with FINITE_CMV (same closure convention).
%   Matrix-valued OPUC: Damanik, Pushnitski and Simon (2008).
%   [59] https://arxiv.org/abs/0711.2703
%   [67] https://doi.org/10.1016/j.laa.2019.11.015
%   Terminal 'pencil' mode returns struct P with P.A, P.B, P.U.
%   Standard [U,PROPERTIES] signature preserves Anymatrix discovery.
properties={'square','unitary','normal','sparse','banded','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(q),q=2;end
validateattributes(q,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3 || isempty(alpha),alpha=zeros(q,q,max(0,n-1));end
if nargin<4 || isempty(beta),beta=eye(q);end
if ~isnumeric(alpha) || any(~isfinite(alpha(:))) ...
        || ~isequal(size(alpha,1),q) || ~isequal(size(alpha,2),q) ...
        || size(alpha,3)~=n-1
    error('block_cmv_pencil:Alpha','ALPHA must be Q-by-Q-by-(N-1).');
end
validateattributes(beta,{'numeric'},{'2d','finite','size',[q q]});
if norm(beta'*beta-eye(q),'fro')>1e-10*max(1,q)
    error('block_cmv_pencil:Beta','BETA must be unitary.');
end
for k=1:n-1
    if norm(alpha(:,:,k),2)>=1
        error('block_cmv_pencil:Alpha','Each ALPHA block must be strictly contractive.');
    end
end
if nargin<5 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('block_cmv_pencil',n,q);
    A=[];B=[];return
end
if ~any(strcmpi(char(mode),{'matrix','pencil'}))
    error('block_cmv_pencil:Mode','Supported modes: matrix, pencil, spectral.');
end
N=n*q;L=speye(N);M=speye(N);
for k=1:n-1
    a=alpha(:,:,k);
    hL=eye(q)-a'*a;hR=eye(q)-a*a';
    rL=sqrtm((hL+hL')/2);rR=sqrtm((hR+hR')/2);
    theta=[a',rL;rR,-a];
    idx=(k-1)*q+(1:2*q);
    if mod(k,2)==1
        L(idx,idx)=theta;
    else
        M(idx,idx)=theta;
    end
end
idx=(n-1)*q+(1:q);
if mod(n,2)==1,L(idx,idx)=beta';else,M(idx,idx)=beta';end
A=M;B=L';
U=L*M;
if strcmpi(char(mode),'pencil')
    U=struct('A',A,'B',B,'U',U);
end
end
