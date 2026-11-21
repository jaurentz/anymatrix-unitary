function Y=um_flag_connection_expmv(X,flag,K,t,ii,jj,twoBlocks)
%UM_FLAG_CONNECTION_EXPMV Projected Lie bracket without a tangent d-by-d matrix.
%   Compute the real normal homogeneous flag Levi-Civita tangent
%   transport exp(t*S), S=-1/2 P_m ad_K, on moving tangent coordinates.
%   K is the real skew n-by-n horizontal geodesic generator from
%   REAL_FLAG_PARALLEL_TRANSPORT, and ii,jj list its off-block pairs.
%   The algorithm never constructs S, whose order d can be O(n^2).
%
%   For V in C^d, assemble the complex-linear skew-coordinate matrix
%      E(ii,jj)=V,  E(jj,ii)=-V  (NO conjugation),
%   and return S*V = -(K*E-E*K)(ii,jj)/2.
%   Each action takes O(n^3) dense arithmetic and O(n^2)
%   working storage, versus O(d^2) for generic dense S.
%   For two blocks, [m,m] lies in the isotropy algebra: S=0.
%
%   A degree-24 scaled Taylor polynomial approximates exp(t*S):
%      m=max(1,ceil(abs(t)*norm(K,'fro')))
%      P24(hS)=sum_{k=0}^{24}(hS)^k/k!, h=t/m.
%   Since ||S||_2 <= ||K||_2 <= ||K||_F, ||hS||<=1.
%   The *exact arithmetic truncation* error bound is at most
%   m*eps*(1+eps)^(m-1) with eps=exp(1)/factorial(25).
%   This DOES NOT include floating roundoff, and the returned
%   truncated action is NOT mathematically exactly orthogonal.
%   The underlying exp(t*S) is orthogonal in real arithmetic.
%
%   FLAG: size, notransp, transp, ctransp; for real S skew,
%   ordinary transpose and conjugate transpose both reverse time.
%   Complex multi-RHS retain linearity; no polar projection.
%
%   [168] https://doi.org/10.1137/24M1684189
n=size(K,1);
d=numel(ii);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_flag_connection_expmv:Flag','FLAG must be text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[d d];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_flag_connection_expmv:Flag', ...
        'Use notransp, transp, ctransp, or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=d
    error('um_flag_connection_expmv:Dimension', ...
        'Input must have the real flag tangent dimension rows.');
end
if issparse(X),X=full(X);end
if strcmp(flag,'notransp'),tau=t;else,tau=-t;end
Y=X;
if tau==0||twoBlocks||~any(K(:)),return;end
knorm=norm(K,'fro');
m=max(1,ceil(abs(tau)*knorm));
if ~isfinite(m)||m>4096
    error('um_flag_connection_expmv:Steps', ...
        'Taylor step count exceeds 4096.');
end
h=tau/m;
upper=sub2ind([n n],ii,jj);
lower=sub2ind([n n],jj,ii);
for step=1:m
    V=Y;
    Z=Y;
    for k=1:24
        V=(h/k)*apply_connection(V);
        Z=Z+V;
    end
    Y=Z;
end

    function W=apply_connection(V)
        W=zeros(size(V),'like',V);
        for rhs=1:size(V,2)
            E=zeros(n,n,'like',V);
            E(upper)=V(:,rhs);
            E(lower)=-V(:,rhs);
            C=K*E-E*K;
            W(:,rhs)=-.5*C(upper);
        end
    end
end
