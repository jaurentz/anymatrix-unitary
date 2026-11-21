function [nodes,info]=um_gauss_sturm_nodes(diagonal,offdiagonal)
%UM_GAUSS_STURM_NODES Gaussian/Jacobi tridiagonal nodes without eigenvectors.
%   [X,INFO]=UM_GAUSS_STURM_NODES(D,E) computes the ordered eigenvalues
%   of the symmetric Jacobi tridiagonal J=diag(D)+diag(E,1)+diag(E,-1),
%   with STRICTLY positive finite E, by inertia-based Sturm bisection.
%
%   A real symmetric tridiagonal's characteristic-inertia count
%   can be evaluated in O(n) via a scalar LDL' recurrence:
%      p_1=d_1-x, p_i=d_i-x-e_(i-1)^2/p_(i-1),
%      # {eigenvalues <=x} = # of negative pivots (with a defined
%      infinitesimal convention for zero pivots).
%   Bisection brackets the k-th eigenvalue using this count.
%   Endpoints initially use rigorous REAL-ARITHMETIC Gershgorin
%   inclusion disks; numerical signs/counts are still subject to
%   finite-precision rounding and are NOT interval-certified.
%
%   Repeating at most 80 steps per node gives O(n^2*80) work and
%   O(n) storage, with no dense J, eigenvectors, or Gauss Q.
%   This is an O(n^2*log(1/u)) eigenvalue-only approach, NOT
%   a fast Hermite/Laguerre polynomial transform.
%   The tiny zero-pivot perturbation is a numerical guard, not a
%   mathematically rigorous finite-precision eigenvalue enclosure.
%
%   INFO reports iteration counts, approximate left/right brackets,
%   their widths, Gershgorin endpoints, maximum width, finite
%   convergence and explicit no-dense-eigenvector metadata.
%   Classical sources: Sturm oscillation; Golub--Welsch (1969).
d=double(diagonal(:));
e=double(offdiagonal(:));
n=numel(d);
if n<1||numel(e)~=max(0,n-1) ...
        ||any(~isfinite(d))||any(~isfinite(e))||any(e<=0)
    error('um_gauss_sturm_nodes:Jacobi', ...
        'Require finite Jacobi diagonal and positive offdiagonals.');
end
if n==1
    nodes=d;
    info=struct('iterations',0,'left_brackets',d,'right_brackets',d, ...
        'bracket_widths',0,'maximum_bracket_width',0, ...
        'gershgorin_interval',d, 'converged',true, ...
        'dense_eigenvectors_constructed',false, ...
        'algorithm','Sturm tridiagonal bisection, O(n^2 log(1/u))');
    return
end
r=[0;e]+[e;0];
scale=max([1;abs(d);abs(e)]);
pad=16*eps(scale);
left=min(d-r)-pad;
right=max(d+r)+pad;
if ~isfinite(left)||~isfinite(right)||right<=left
    error('um_gauss_sturm_nodes:Bounds', ...
        'Gershgorin interval cannot be represented.');
end
tiny=1e-100*scale;
nodes=zeros(n,1);
lefts=zeros(n,1);
rights=zeros(n,1);
iters=zeros(n,1);
converged=true;
for k=1:n
    a=left;b=right;
    for it=1:80
        mid=a+(b-a)/2;
        if mid==a||mid==b
            break
        end
        count=sturmcount(d,e,mid,tiny);
        if count<k
            a=mid;
        else
            b=mid;
        end
        if b-a<=4*eps(max([1,abs(a),abs(b)]))
            break
        end
    end
    nodes(k)=a+(b-a)/2;
    lefts(k)=a;
    rights(k)=b;
    iters(k)=it;
    if b-a>32*eps(max([1,abs(a),abs(b)]))
        converged=false;
    end
end
if ~converged||any(diff(nodes)<=0)||any(~isfinite(nodes))
    error('um_gauss_sturm_nodes:Resolution', ...
        'Sturm bisection failed to resolve distinct Jacobi nodes.');
end
widths=rights-lefts;
info=struct('iterations',iters,'left_brackets',lefts, ...
    'right_brackets',rights,'bracket_widths',widths, ...
    'maximum_bracket_width',max(widths), ...
    'gershgorin_interval',[left right], ...
    'converged',converged, ...
    'dense_eigenvectors_constructed',false, ...
    'algorithm','Sturm tridiagonal bisection, O(n^2 log(1/u))', ...
    'notes',['Finite precision inertia/bisection brackets are ', ...
        'numerical estimates, not certified interval eigenvalue ', ...
        'enclosures; no dense orthogonal Gauss matrix was formed.']);
end

function count=sturmcount(d,e,x,tiny)
pivot=d(1)-x;
if abs(pivot)<tiny,pivot=-tiny;end
count=double(pivot<0);
for j=2:numel(d)
    pivot=d(j)-x-(e(j-1)^2)/pivot;
    if abs(pivot)<tiny,pivot=-tiny;end
    count=count+double(pivot<0);
end
end
