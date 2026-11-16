function [src,dst,reverse,degree,starts]=um_graph_arc_data(A)
%UM_GRAPH_ARC_DATA Canonical directed-arc ordering of a simple graph.
%   Arcs (v,w) are sorted first by v and then by w. REVERSE maps every
%   directed arc to (w,v). Reject isolated vertices and loops.
validateattributes(A,{'numeric'},{'2d','real','finite','nonnegative'});
n=size(A,1);
if size(A,2)~=n || n<2 || any(diag(A)~=0) ...
        || norm(A-A','fro')>1e-12*max(1,norm(A,'fro'))
    error('um_graph_arc_data:Adjacency', ...
          'Need a symmetric real nonnegative zero-diagonal matrix (n>=2).');
end
[i,j]=find(triu(A,1));
src=[i;j];dst=[j;i];
[~,p]=sortrows([src,dst],[1,2]);src=src(p);dst=dst(p);
degree=accumarray(src,1,[n 1]);
if any(degree==0)
    error('um_graph_arc_data:Isolated','Every vertex must have an edge.');
end
m=numel(src);
starts=[1;cumsum(degree)+1];
index=sparse(src,dst,1:m,n,n);
reverse=full(index(sub2ind([n n],dst,src)));
reverse=reverse(:);
end
