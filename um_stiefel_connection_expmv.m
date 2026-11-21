function Y=um_stiefel_connection_expmv(X,flag,n,p,A,B,t,ii,jj)
%UM_STIEFEL_CONNECTION_EXPMV Matrix-free canonical Stiefel transport action.
%   Approximates exp(t*S)*X in moving orthonormal tangent coordinates,
%   S=-P_horizontal[ K, . ]/2, K=[A,-B';B,0].
%   Neither S (d-by-d) nor K (n-by-n) is formed. An individual
%   commutator action costs O(n*p^2) for p<<n, using a p-by-p skew
%   block and an (n-p)-by-p block. Storage O(n*p+numel(X)).
%
%   Every action is a scaled Taylor matrix exponential (degree 24):
%   split t into m=max(1,ceil(abs(t)*||K||_F)) steps. Since
%   ||S||_2 <= ||K||_2 <= ||K||_F, each step has ||h*S||<=1.
%   In exact arithmetic the one-step Taylor truncation error is
%   <= exp(1)/25!, and the m-step error is bounded by
%   m*eps24*(1+eps24)^(m-1). For m<=4096 this is <1e-20.
%   Floating roundoff is NOT covered by this exact-arithmetic
%   bound: this action is numerically near-isometric, NOT an
%   exactly orthogonal finite polynomial. The original dense
%   expm(t*S) is the exact mathematical reference.
%
%   FLAG 'notransp','ctransp','transp','size'. Since S is real
%   skew-symmetric, the transpose/adjoint actions use -t.
%   Complex multiple RHS use complex-linear (nonconjugated)
%   horizontal skew-form coordinates. [168]
%   https://doi.org/10.1137/24M1684189
d=numel(ii);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_stiefel_connection_expmv:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[d d];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_stiefel_connection_expmv:Flag', ...
        'Use notransp, transp, ctransp, or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=d
    error('um_stiefel_connection_expmv:Dimension', ...
        'X must have N*P-P*(P+1)/2 rows.');
end
if issparse(X),X=full(X);end
if strcmp(flag,'notransp'),tau=t;else,tau=-t;end
knorm=hypot(norm(A,'fro'),sqrt(2)*norm(B,'fro'));
m=max(1,ceil(abs(tau)*knorm));
if ~isfinite(m)||m>4096
    error('um_stiefel_connection_expmv:Steps', ...
        'Scaled Taylor bound requires at most 4096 time steps.');
end
Y=X;
if tau==0||knorm==0||p==1,return;end
h=tau/m;
top=jj<=p;
cross=~top;
indUp=sub2ind([p p],ii(top),jj(top));
indDown=sub2ind([p p],jj(top),ii(top));
if n>p
    indMixed=sub2ind([n-p p],jj(cross)-p,ii(cross));
else
    indMixed=[];
end
for step=1:m
    Z=Y;R=Y;
    for k=1:24
        R=(h/k)*apply_connection(R);
        Z=Z+R;
    end
    Y=Z;
end

    function W=apply_connection(V)
        W=zeros(size(V),'like',V);
        for rhs=1:size(V,2)
            E11=zeros(p,p,'like',V);
            E11(indUp)=V(top,rhs);
            E11(indDown)=-V(top,rhs);
            E21=zeros(n-p,p,'like',V);
            if n>p,E21(indMixed)=-V(cross,rhs);end
            W11=A*E11-E11*A-B'*E21+E21'*B;
            W21=B*E11-E21*A;
            W(top,rhs)=-.5*W11(indUp);
            if n>p,W(cross,rhs)=.5*W21(indMixed);end
        end
    end
end
