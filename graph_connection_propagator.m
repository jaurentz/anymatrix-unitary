function [U,properties]=graph_connection_propagator(n,q,edges,weights,links,t,output,cayley_steps)
%GRAPH_CONNECTION_PROPAGATOR Unitary graph connection-Laplacian flow.
%   U=GRAPH_CONNECTION_PROPAGATOR(N,Q,EDGES,WEIGHTS,LINKS,T,OUTPUT)
%   returns exp(-1i*T*L) of a Hermitian positive-semidefinite
%   graph connection Laplacian L, order N*Q.
%   Each undirected edge e is oriented i=EDGES(e,1),j=EDGES(e,2);
%   LINKS(:,:,e) is a Q-by-Q unitary transport mapping fiber j to i.
%   L_ii += w_e I, L_jj += w_e I, L_ij -= w_e R_ij,
%   L_ji -= w_e R_ij'. Edges are distinct, without self loops,
%   weights nonnegative, and links unitary. Consequently
%      x'*L*x = sum_e w_e*norm(x_i-R_ij*x_j)^2 >= 0.
%   Under local unitary frames G_i, replacing R_ij with
%   G_i*R_ij*G_j' conjugates L and exp(-itL) by diag(G_i).
%
%   By default the graph is a unit-weight N-cycle (N>=3) with
%   Q independent holonomy phases .33*(1:Q), concentrated on
%   the closing edge (N,1). For channel a the exact L spectrum is
%      lambda_(a,k)=2-2*cos((2*pi*k+.33*a)/N), k=0,...,N-1.
%   Thus the default output='spectral' gives exact finite phase
%   multiset exp(-1i*T*lambda_(a,k)) and trace moments.
%   For arbitrary supplied graphs the generic analytic spectral
%   descriptor is UNKNOWN; eig(L) is not silently called analytic.
%
%   OUTPUT='matrix' (default) constructs TRUE expm(-1i*T*L),
%   limited to moderate order <=256. OUTPUT='laplacian' returns
%   the sparse Hermitian generator and its edge/frame metadata.
%   OUTPUT='cayley_operator' returns a LABELED approximation struct
%   with APPLY(X,FLAG), dimension and CAYLEY_STEPS. Its operator is
%      [ (I+itL/(2s))\(I-itL/(2s)) ]^s
%   rather than exp(-itL). It is exactly unitary in exact arithmetic,
%   uses cached sparse factorizations and has O(s^-2) fixed-t error.
%   It is NOT advertised as exact matrix-free exponential action.
%   A generic sparse Krylov expmv likewise does not preserve exact
%   unitarity, so it is not mislabeled as a strict unitary.
%
%   [123] Singer and Wu (2012):
%   https://doi.org/10.1002/cpa.21395
%   [98] Crane, Desbrun, Schroder (2010):
%   https://doi.org/10.1111/j.1467-8659.2010.01761.x
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'finite'});
if nargin<2||isempty(q),q=2;end
validateattributes(q,{'numeric'},{'scalar','integer','positive','finite'});
defaultModel=(nargin<3||isempty(edges)) ...
    &&(nargin<4||isempty(weights)) ...
    &&(nargin<5||isempty(links));
if nargin<3||isempty(edges)
    edges=[(1:n-1)',(2:n)';n,1];
end
if ~isnumeric(edges)||~isreal(edges)||size(edges,2)~=2 ...
        ||any(~isfinite(edges(:)))||any(edges(:)<1) ...
        ||any(edges(:)>n)||any(edges(:)~=round(edges(:))) ...
        ||any(edges(:,1)==edges(:,2))
    error('graph_connection_propagator:Edges', ...
        'EDGES must be distinct oriented pairs of different vertices.');
end
edges=double(edges);
m=size(edges,1);
if size(unique(sort(edges,2),'rows'),1)~=m
    error('graph_connection_propagator:Edges', ...
        'Repeated undirected edges are not supported.');
end
if nargin<4||isempty(weights),weights=ones(m,1);end
if ~isnumeric(weights)||~isreal(weights)||~isvector(weights) ...
        ||numel(weights)~=m||any(~isfinite(weights(:))) ...
        ||any(weights(:)<0)
    error('graph_connection_propagator:Weights', ...
        'WEIGHTS must contain M finite nonnegative edge weights.');
end
weights=weights(:);
if nargin<5||isempty(links)
    links=repmat(eye(q),1,1,m);
    if defaultModel
        links(:,:,m)=diag(exp(1i*.33*(1:q)));
    end
end
if ~isnumeric(links)||any(~isfinite(links(:))) ...
        ||size(links,1)~=q||size(links,2)~=q ...
        ||size(links,3)~=m
    error('graph_connection_propagator:Links', ...
        'LINKS must be finite Q-by-Q-by-M transport blocks.');
end
for e=1:m
    R=links(:,:,e);
    if norm(R'*R-eye(q),'fro')>2e-10*max(1,q)
        error('graph_connection_propagator:Links', ...
            'Every link must be unitary.');
    end
end
if nargin<6||isempty(t),t=1;end
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<7||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('graph_connection_propagator:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if nargin<8||isempty(cayley_steps),cayley_steps=16;end
validateattributes(cayley_steps,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
N=n*q;
if N>flintmax||~isfinite(N)
    error('graph_connection_propagator:Order','Dimension exceeds index range.');
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('graph_connection_propagator',n);
    U.dimension=N;
    U.reference='Singer and Wu (2012), vector diffusion maps and the connection Laplacian';
    if defaultModel
        phi=.33*(1:q)';
        k=(0:n-1);
        lambda=2-2*cos((2*pi*k+phi)/n);
        z=exp(-1i*t*lambda(:));
        U.finite.kind='exact';
        U.finite.eigenvalues=@()z;
        U.finite.moments=@(orders)um_graph_connection_moments(z,orders);
        U.finite.atoms=struct('angles',angle(z),'weights',ones(N,1)/N);
        U.notes=['Exact twisted-cycle fiber holonomy spectral phases; ' ...
            'the trace measure is uniform over all N*Q eigenphases.'];
    else
        U.notes=['Generic weighted graph connection Laplacian: ' ...
            'no asserted closed-form finite trace eigenphase law.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','laplacian','cayley_operator'}))
    error('graph_connection_propagator:Mode', ...
        'Choose matrix, laplacian, cayley_operator, or spectral.');
end
L=spalloc(N,N,4*m*q*q);
for e=1:m
    i=edges(e,1);j=edges(e,2);w=weights(e);
    I=(i-1)*q+(1:q);J=(j-1)*q+(1:q);
    R=links(:,:,e);
    L(I,I)=L(I,I)+w*speye(q);
    L(J,J)=L(J,J)+w*speye(q);
    L(I,J)=L(I,J)-w*R;
    L(J,I)=L(J,I)-w*R';
end
switch output
    case 'laplacian'
        U=struct('L',L,'dimension',N,'fiber_dimension',q, ...
            'edges',edges,'weights',weights,'links',links, ...
            'orientation','R_ij maps fiber j to vertex i', ...
            'gauge_convention','Rnew_ij=Gi*Rij*Gj''');
    case 'cayley_operator'
        h=t/cayley_steps;
        Ip=speye(N)+.5i*h*L;
        Im=speye(N)-.5i*h*L;
        plusFactor=decomposition(Ip,'lu');
        minusFactor=decomposition(Im,'lu');
        F=@(X,flag)um_graph_connection_cayley(X,flag,N,cayley_steps, ...
            plusFactor,minusFactor);
        U=struct('apply',F,'dimension',N,'steps',cayley_steps, ...
            'approximation','unitary_Cayley_rational_not_exact_expm', ...
            'notes',['Sparse Cayley/implicit-midpoint rational approximate ' ...
                'propagator with O(steps^-2) fixed-time error; not the ' ...
                'true expm(-itL) matrix.']);
    otherwise
        if N>256
            error('graph_connection_propagator:Order', ...
                'Exact dense expm output limited to N*Q<=256.');
        end
        U=expm(-1i*t*full(L));
end
end
