function Y=um_wigner_skew_expmv(X,flag,ell,alpha,beta,gamma)
%UM_WIGNER_SKEW_EXPMV Matrix-free Wigner D action for arbitrary Euler angles.
%   Exact mathematical matrix is
%        D=diag(exp(-1i*alpha*m))*exp(beta*K)*diag(exp(-1i*gamma*m)),
%   with m=-ell:ell, n=2*ell+1 and real skew-tridiagonal
%        K(j,j+1)=sqrt((ell-m_j)*(ell+m_j+1))/2,
%        K(j+1,j)=-K(j,j+1).
%
%   A degree-32 Taylor action applied in ceil(|betaReduced|*ell/2)
%   substeps uses only the O(n) tridiagonal generator. In the spin
%   ell representation, ||K||_2=ell; each substep has ||h*K||_2<=2.
%   In exact arithmetic the Taylor remainder per substep is bounded
%   by epsilon=exp(2)*2^33/33!, and the m-step global absolute
%   operator bound is <=m*epsilon*(1+epsilon)^(m-1).
%   (This is a mathematical truncation bound, NOT a floating-point
%   error certification. Roundoff accumulates in repeated actions.)
%   The SU(2) action is 4*pi periodic in beta, also valid for integer
%   spins. Reducing beta modulo 4*pi avoids unnecessary long stepping.
%
%   This is a rigorously bounded *approximate* exponential action:
%   NO dense n-by-n D or K and NO exact finite unitarity of the
%   computed floating-point action are claimed. Complexity is
%   O(n*32*ceil(|betaReduced|*ell/2)) per RHS, O(n) scratch per
%   RHS; at fixed reduced beta, O(ell^2) per irreducible block,
%   not a subquadratic single-degree FFT/butterfly Wigner transform.
%   Supports complex-linear forward, TRANSPOSE and ADJOINT multi-RHS.
%
%   Sources: angular-momentum tridiagonal generator and
%   exp(-i beta Jy) Wigner convention, as in wigner_rotation.m.
n=round(2*ell)+1;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_wigner_skew_expmv:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_wigner_skew_expmv:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_wigner_skew_expmv:Dimension', ...
        'Input must have 2*ell+1 rows.');
end
m=(-ell:ell)';
pa=exp(-1i*alpha*m);
pg=exp(-1i*gamma*m);
switch flag
    case 'notransp'
        Y=bsxfun(@times,pg,full(X));
        signBeta=1;
    case 'transp'
        Y=bsxfun(@times,pa,full(X));
        signBeta=-1;
    otherwise % conjugate transpose
        Y=bsxfun(@times,conj(pa),full(X));
        signBeta=-1;
end
if ell>0
    b=mod(beta+2*pi,4*pi)-2*pi;
    steps=max(1,ceil(abs(b)*ell/2));
    if steps>8192
        error('um_wigner_skew_expmv:Steps', ...
            'Taylor action requires more than 8192 steps.');
    end
    h=(signBeta*b)/steps;
    c=.5*sqrt((ell-m(1:end-1)).*(ell+m(1:end-1)+1));
    for step=1:steps
        T=Y;
        for k=1:32
            KT=zeros(size(T),'like',T);
            KT(1:end-1,:)=bsxfun(@times,c,T(2:end,:));
            KT(2:end,:)=KT(2:end,:)-bsxfun(@times,c,T(1:end-1,:));
            T=(h/k)*KT;
            Y=Y+T;
        end
    end
end
switch flag
    case 'notransp'
        Y=bsxfun(@times,pa,Y);
    case 'transp'
        Y=bsxfun(@times,pg,Y);
    otherwise
        Y=bsxfun(@times,conj(pg),Y);
end
end
