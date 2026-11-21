function [G,info]=um_wigner_jacobi(L,m,n,x)
%UM_WIGNER_JACOBI L2(dx)-orthonormal fixed-index Wigner d latitude modes.
%   G=UM_WIGNER_JACOBI(L,M,N,X) returns the Lx(L-j0) (if numel(X)=L)
%   real sampled angular functions G_(ell)(x), ell=j0,...,L-1,
%   j0=max(abs(m),abs(n)), for two INTEGER Wigner indices m,n:
%       G_ell(x)=sqrt((2ell+1)/2)*d^ell_(m,n)(acos(x)).
%   For the repo's K=(diag(c,1)-diag(c,-1))/2, with
%   d^ell=expm(theta*K), this is EXACTLY
%       sign * (1-x)^(a/2)*(1+x)^(b/2)*p_(ell-j0)(x)/sqrt(mu0),
%       a=abs(m-n), b=abs(m+n), sign=(-1)^max(0,m-n),
%   where p_k is positive-leading Jacobi orthonormal for the
%   NORMALIZED measure ((1-x)^a(1+x)^b/mu0)dx and
%       mu0=2^(a+b+1)*B(a+1,b+1).
%   Thus integral_-1^1 G_ell*G_ell' dx=delta_(ell,ell').
%   This identity was independently verified for each fixed pair
%   -5<=m,n<=5 against WIGNER_ROTATION at four arbitrary x values.
%
%   The standard spin harmonics in our convention are
%       {}_sY_(ell,m)(theta,phi)=
%         (-1)^s*G_ell(cos(theta))*exp(i*m*phi)/sqrt(2*pi),
%   with n=-s. Peter-Weyl functions satisfy
%       sqrt(2ell+1)*D^ell_(m,n)(a,theta,g)=
%         sqrt(2)*G_ell(cos(theta))*exp(-i*m*a-i*n*g).
%   This helper computes neither a unitary square transform
%   nor a complete angular sampling isometry by itself.
%
%   Positive-leading normalized Jacobi recurrence is evaluated
%   directly in O(numel(x)*(L-j0)), with no dense expm.
%   Restricted L<=64 (no arbitrary-large-degree stability claim).
%
%   [165] Huffenberger-Wandelt, ApJS 189 (2010):
%   https://doi.org/10.1088/0067-0049/189/2/255
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',64});
validateattributes(m,{'numeric'}, ...
    {'scalar','integer','real','finite'});
validateattributes(n,{'numeric'}, ...
    {'scalar','integer','real','finite'});
if max(abs(m),abs(n))>=L
    error('um_wigner_jacobi:Order','Both fixed indices must have absolute value < L.');
end
if ~isnumeric(x)||~isreal(x)||~isvector(x)|| ...
        isempty(x)||any(~isfinite(x(:)))||any(abs(x(:))>=1)
    error('um_wigner_jacobi:Nodes','X must be finite in the open interval (-1,1).');
end
x=x(:);j0=max(abs(m),abs(n));N=L-j0;
a=abs(m-n);b=abs(m+n);
sgn=(-1)^max(0,m-n);
logmu=(a+b+1)*log(2)+gammaln(a+1)+ ...
    gammaln(b+1)-gammaln(a+b+2);
fac=sgn*exp(.5*a*log1p(-x)+.5*b*log1p(x)-.5*logmu);
G=zeros(numel(x),N);
p0=ones(numel(x),1);pm=0*p0;
G(:,1)=fac.*p0;
if N>1
    diagJac=zeros(N,1);offJac=zeros(N-1,1);
    diagJac(1)=(b-a)/(a+b+2);
    k=(1:N-1)';
    diagJac(2:end)=(b^2-a^2)./ ...
        ((2*k+a+b).*(2*k+a+b+2));
    offJac(1)=2/(a+b+2)*sqrt((1+a)*(1+b)/(a+b+3));
    if N>2
        k=(2:N-1)';
        offJac(2:end)=2./(2*k+a+b).* ...
            sqrt(k.*(k+a).*(k+b).*(k+a+b)./ ...
            ((2*k+a+b-1).*(2*k+a+b+1)));
    end
    for k=1:N-1
        if k==1, old=0;else,old=offJac(k-1)*pm;end
        pnext=((x-diagJac(k)).*p0-old)/offJac(k);
        G(:,k+1)=fac.*pnext;
        pm=p0;p0=pnext;
    end
end
if nargout>1
    info=struct('m',m,'n',n,'ell_min',j0, ...
        'ell_values',(j0:L-1)', ...
        'jacobi_alpha',a,'jacobi_beta',b,'phase',sgn, ...
        'log_jacobi_mass',logmu, ...
        'normalization','G_ell=sqrt((2ell+1)/2)*Wigner-d_(m,n)');
end
end
