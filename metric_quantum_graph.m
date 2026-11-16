function [U,properties]=metric_quantum_graph(n,lengths,k,coins,mode)
%METRIC_QUANTUM_GRAPH Unitary bond-scattering map on a metric graph.
%   U=METRIC_QUANTUM_GRAPH(N,LENGTHS,K,COINS) returns the unitary
%       U(K)=D(K)*S, D(K)=diag(exp(i*K*length_b)),
%   on the 2m directed arcs of an undirected simple graph.
%   LENGTHS is N-by-N symmetric, zero-diagonal, strictly positive on
%   edges, in physical length units; zero entries denote missing edges.
%   K is a real wavenumber. Default metric graph is an irregular cycle.
%   S=C*R, R reverses arcs and C applies vertex-local scattering coins.
%   COINS='kirchhoff' (default) uses Neumann-Kirchhoff vertex scatterers
%       C_v=(2/d_v)*ones(d_v)-eye(d_v).
%   COINS='fourier', 'identity', or a cell array of arbitrary unitary
%   local vertex matrices is also accepted. Arc neighbors sorted ascending.
%
%   MODE='operator' returns a matrix-free fixed-K unitary application.
%   MODE='secular' returns a struct with U(z), T(z)=I-U(z), det(T(z))
%   as analytic matrix functions of complex z, including nonreal z where
%   U(z) is generally not unitary. The real secular wavenumbers satisfy
%       det(I-U(K))=0.
%   For incommensurate metric lengths this is a nonlinear, generally
%   transcendental eigenvalue problem, not one finite unitary EVP.
%   Unlike UNITARY_DELAY_PENCIL there is no generic finite polynomial lift
%   for real incommensurate lengths.
%   MODE='spectral' describes eigenphases of the fixed real-K matrix only.
%
%   Kottos & Smilansky, Ann. Phys. 274 (1999), 76-124.
%   https://doi.org/10.1006/aphy.1999.5904
%   Gutkin & Smilansky, J. Phys. A 34 (2001), 6061-6068.
%   https://doi.org/10.1088/0305-4470/34/31/301
properties={'square','unitary','normal','complex','sparse','scalable', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',2});
if nargin<2 || isempty(lengths)
    d=1+sqrt(2)*(1:n-1)'/(n+1);
    lengths=sparse(1:n-1,2:n,d,n,n);
    if n>2,lengths(1,n)=1+sqrt(3)/5;end
    lengths=lengths+lengths';
end
validateattributes(lengths,{'numeric'},{'size',[n n]});
[src,dst,reverse,~,starts]=um_graph_arc_data(lengths);
N=numel(src);
if nargin<3 || isempty(k),k=1;end
validateattributes(k,{'numeric'},{'scalar','real','finite'});
if nargin<4 || isempty(coins),coins='kirchhoff';end
if nargin<5 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('metric_quantum_graph:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',N);
    U.family='metric_quantum_graph';
    U.reference='Kottos and Smilansky, Ann. Phys. 274 (1999)';
    U.notes=['Eigenphases of U(k) at a fixed real k are not the ' ...
        'quantum-graph wavenumbers solving det(I-U(k))=0. ' ...
        'Real incommensurate lengths require a transcendental NEP.'];
    return
end
if ~any(strcmp(mode,{'matrix','operator','secular'}))
    error('metric_quantum_graph:Mode', ...
          'Supported modes: matrix, operator, secular, spectral.');
end
if (ischar(coins)||isstring(coins)) && strcmpi(char(coins),'kirchhoff')
    coins='grover';
end
C=um_graph_coin_matrix(lengths,src,dst,starts,coins,false);
L=full(lengths(sub2ind([n n],src,dst)));
if strcmp(mode,'secular')
    reverse=reverse(:);
    R=sparse(reverse,1:N,1,N,N);
    base=C*R;
    F=@(z)spdiags(exp(1i*z*L),0,N,N)*base;
    U=struct('U',F,'T',@(z)speye(N)-F(z), ...
             'det',@(z)det(full(speye(N)-F(z))), ...
             'dimension',N,'lengths',L, ...
             'notes','U(z) is unitary only for real z; T(z) is an analytic NEP.');
    return
end
phase=exp(1i*k*L);
if strcmp(mode,'operator')
    U=@(X,flag)um_graph_walk_apply(X,flag,C,reverse,phase,'metric');
    return
end
R=sparse(reverse,1:N,1,N,N);
U=spdiags(phase,0,N,N)*(C*R);
end
