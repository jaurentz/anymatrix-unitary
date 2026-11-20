function [U,properties]=block_cmv_pencil(n,q,alpha,beta,mode)
%BLOCK_CMV_PENCIL Finite unitary block-CMV matrix and strict pencil.
%   U=BLOCK_CMV_PENCIL(N,Q,ALPHA,BETA) is an N*Q-by-N*Q
%   block-CMV UNITARY MATRIX with noncommuting matrix-valued
%   Verblunsky coefficients ALPHA(:,:,j), j=1,...,N-1.
%   Every Q-by-Q alpha_j is strictly contractive and the terminal
%   Q-by-Q BETA is unitary. The alternating core factors L,M satisfy
%       U=L*M,   Theta_j=[alpha_j',rho_L;rho_R,-alpha_j],
%       rho_L=sqrt(I-alpha_j'*alpha_j),
%       rho_R=sqrt(I-alpha_j*alpha_j').
%   The distinct left and right defect matrices must not be swapped.
%   Mode 'pencil' returns strict factors A=M, B=L', U=L*M;
%   thus B\A=U (the matrix output is already genuinely unitary).
%   Mode 'operator' returns O(X,FLAG) for matrix-free forward/adjoint
%   O(N*Q^2) actions via local block cores and boundary mixing.
%   Mode 'quadrature' computes NUMERICAL N*Q eigenphases and the
%   Q-by-Q matrix-valued cyclic spectral weights
%       W_j=V(1:Q,j)*V(1:Q,j)'.
%   Sum W_j=I_Q; moments sum lambda_j^k*W_j are the upper-left Q-by-Q
%   block of U^k. This is NOT the uniform scalar trace measure,
%   which has weights 1/(N*Q). Repeated-eigenvalue projector
%   weights are invariant only after aggregating coincident nodes.
%   'spectral' retains the UNKNOWN analytic generic phase descriptor.
%   Q=1 reduces to FINITE_CMV with the same terminal convention.
%   [59] https://arxiv.org/abs/0711.2703
%   [67] https://doi.org/10.1016/j.laa.2019.11.015
properties={'square','unitary','normal','sparse','banded','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(q),q=2;end
validateattributes(q,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3||isempty(alpha),alpha=zeros(q,q,max(0,n-1));end
if nargin<4||isempty(beta),beta=eye(q);end
if ~isnumeric(alpha)||any(~isfinite(alpha(:))) ...
        ||size(alpha,1)~=q||size(alpha,2)~=q ...
        ||size(alpha,3)~=n-1
    error('block_cmv_pencil:Alpha','ALPHA must be Q-by-Q-by-(N-1).');
end
validateattributes(beta,{'numeric'},{'2d','finite','size',[q q]});
if norm(beta'*beta-eye(q),'fro')>1e-10*max(1,q)
    error('block_cmv_pencil:Beta','BETA must be unitary.');
end
for k=1:n-1
    if norm(alpha(:,:,k),2)>=1
        error('block_cmv_pencil:Alpha','Each ALPHA must be strictly contractive.');
    end
end
if nargin<5||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('block_cmv_pencil:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    U=ep_spectral_reference('block_cmv_pencil',n,q);
    return
end
if ~any(strcmp(mode,{'matrix','pencil','operator','implicit','quadrature'}))
    error('block_cmv_pencil:Mode', ...
        'Supported modes: matrix, pencil, operator, quadrature, spectral.');
end
N=n*q;
cores=cell(n-1,1);
for k=1:n-1
    a=alpha(:,:,k);
    hL=eye(q)-a'*a;
    hR=eye(q)-a*a';
    rL=sqrtm((hL+hL')/2);
    rR=sqrtm((hR+hR')/2);
    cores{k}=[a',rL;rR,-a];
end
if any(strcmp(mode,{'operator','implicit'}))
    U=@(X,flag)um_block_cmv_apply(X,flag,n,q,cores,beta);
    return
end
L=speye(N);M=speye(N);
for k=1:n-1
    idx=(k-1)*q+(1:2*q);
    if mod(k,2)==1
        L(idx,idx)=cores{k};
    else
        M(idx,idx)=cores{k};
    end
end
idx=(n-1)*q+(1:q);
if mod(n,2)==1,L(idx,idx)=beta';else,M(idx,idx)=beta';end
U=L*M;
if strcmp(mode,'pencil')
    U=struct('A',M,'B',L','U',U);
elseif strcmp(mode,'quadrature')
    if N>128
        error('block_cmv_pencil:QuadratureSize', ...
            'Dense numerical spectral projectors limited to N*Q<=128.');
    end
    [V,D]=eig(full(U));
    nodes=diag(D);
    [~,order]=sort(angle(nodes),'ascend');
    nodes=nodes(order);
    V=V(:,order);
    W=complex(zeros(q,q,N));
    for j=1:N
        v=V(1:q,j);
        W(:,:,j)=v*v';
    end
    U=struct('kind','numerical','dimension',N,'fiber_dimension',q, ...
        'nodes',nodes,'weights',W, ...
        'trace_weights',ones(N,1)/N, ...
        'cyclic_moments',@(k)um_block_cmv_moments(nodes,W,k,'cyclic'), ...
        'trace_moments',@(k)um_block_cmv_moments(nodes,W,k,'trace'), ...
        'notes',['Numerical block cyclic weights from eig(U); ' ...
            'not an analytic finite phase reference. ' ...
            'For repeated nodes aggregate spectral projectors.']);
end
end
