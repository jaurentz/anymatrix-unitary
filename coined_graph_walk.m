function [U,properties]=coined_graph_walk(n,adjacency,coins,mode)
%COINED_GRAPH_WALK Flip-flop coined quantum walk on an arbitrary graph.
%   U=COINED_GRAPH_WALK(N,ADJACENCY,COINS) constructs the 2m-by-2m
%   unitary U=S*C on ordered directed arcs (vertex,neighbor), where m
%   is the number of edges of an undirected graph. ADJACENCY is an
%   N-by-N symmetric nonnegative zero-diagonal matrix with no isolates;
%   positive values can be treated as edge weights.
%
%   Local unitary coins can be: 'grover' (default, weighted by edge
%   weights), 'fourier', 'identity', or a cell array of N local d_v-by-d_v
%   unitary matrices, with arc neighbors sorted in ascending order.
%   S reverses each directed edge, hence S=S'=S^(-1).
%   A=COINED_GRAPH_WALK(...,'pencil') returns A.A=C, A.B=S',
%   A.U=S*C, a strict generalized unitary pencil.
%   MODE='operator' gives an O(sum d_v^2) sparse-coin matvec handle.
%   MODE='spectral' returns unknown finite trace metadata.
%   A graph Fourier basis is not the eigenbasis of this walk.
%
%   D. Emms, S. Severini, R. C. Wilson, E. R. Hancock,
%   Coined quantum walks lift the cospectrality of graphs and trees,
%   Pattern Recognition 42 (2009), 1988-2002.
%   https://doi.org/10.1016/j.patcog.2008.10.025
properties={'square','unitary','normal','sparse','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',2});
if nargin<2 || isempty(adjacency)
    adjacency=sparse(1:n-1,2:n,1,n,n);
    adjacency=adjacency+adjacency';
end
validateattributes(adjacency,{'numeric'},{'size',[n n]});
[src,dst,reverse,~,starts]=um_graph_arc_data(adjacency);
N=numel(src);
if nargin<3 || isempty(coins),coins='grover';end
if nargin<4 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('coined_graph_walk:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',N);
    U.family='coined_graph_walk';
    U.reference='Emms et al., Pattern Recognition 42 (2009)';
    U.notes=['Flip-flop walk with vertex-local unitary coins. ' ...
        'Grover coins admit a normalized-adjacency spectral mapping, ' ...
        'but general local coins have no universal analytic trace law.'];
    return
end
if ~any(strcmp(mode,{'matrix','pencil','operator'}))
    error('coined_graph_walk:Mode', ...
          'Supported modes: matrix, pencil, operator, spectral.');
end
C=um_graph_coin_matrix(adjacency,src,dst,starts,coins,true);
if strcmp(mode,'operator')
    U=@(X,flag)um_graph_walk_apply(X,flag,C,reverse,[], 'coined');
    return
end
S=sparse(reverse,1:N,1,N,N);
U=S*C;
if strcmp(mode,'pencil')
    U=struct('A',C,'B',S','U',U);
end
end
