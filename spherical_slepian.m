function [U,properties]=spherical_slepian(L,cap,output)
%SPHERICAL_SLEPIAN Real orthogonal polar-cap Slepian eigenbasis on S^2.
%   U=SPHERICAL_SLEPIAN(L,CAP) returns the L^2-by-L^2 ANALYSIS
%   matrix U=Q', where columns of Q are orthonormal real eigenvectors
%   of the bandlimited spherical concentration matrix
%       K_ab = integral_{0<=theta<=CAP} Y_a Y_b dOmega.
%   Degrees are ell=0,...,L-1 (L is the BANDLIMIT, not matrix order).
%   CAP is the polar-cap colatitude in [0,pi], default pi/3.
%   The physical basis ordering is: m=0, ell=0,...,L-1, then for
%   m=1,...,L-1, first cos(m*phi) then sin(m*phi), ell=m,...,L-1.
%   Each real harmonic is L2(S^2)-normalized with Condon--Shortley
%   associated Legendre convention. Q synthesis, U=Q' analysis.
%
%   Within each fixed-m harmonic block, an L-point Gauss--Legendre
%   rule on x=cos(CAP),...,1 EXACTLY integrates products P_l^m P_k^m:
%   their product is a polynomial of degree <=2L-2. Cap symmetry
%   makes all cross-m blocks exactly zero. Its concentration factors
%   lie in [0,1], but are NOT eigenvalues/eigenphases of U.
%   Eigenvector signs have largest-magnitude coefficient positive;
%   all eigenvectors are sorted by decreasing concentration (stable
%   original block-index tie break). The repeated cosine/sine
%   eigenvalues use separate block bases, avoiding arbitrary mixing.
%   Full/empty caps use the canonical Q=I.
%
%   OUTPUT 'matrix'/'analysis' (default): Q' (orthogonal);
%          'basis': Q (orthogonal synthesis);
%          'concentration': K (NOT unitary);
%          'parameters': struct with K,Q,eigenvalues, basis labels;
%          'spectral': unknown eigenphases of Q', with concentration
%                      factors explicitly identified as a DIFFERENT law.
%   No generic polygon/mask region or fast structured operator is
%   claimed: those require separate accuracy and complexity studies.
%
%   [156] Simons, Dahlen, Wieczorek, SIAM Rev. 48 (2006), 504--536.
%   https://doi.org/10.1137/S0036144504445765
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'},{'scalar','integer','positive','finite','<=',24});
if nargin<2||isempty(cap),cap=pi/3;end
validateattributes(cap,{'numeric'},{'scalar','real','finite','>=',0,'<=',pi});
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('spherical_slepian:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','analysis','basis','concentration', ...
        'parameters','spectral'}))
    error('spherical_slepian:Mode','Unknown output mode.');
end
n=L^2;
if strcmp(output,'spectral')
    U=unitary_spectral_reference('spherical_slepian',n);
    U.bandlimit=L;U.cap=cap;
    U.reference='Simons--Dahlen--Wieczorek, SIAM Review 48 (2006)';
    U.notes=['Orthogonal Slepian analysis eigenbasis; its eigenphases ' ...
        'are generally unknown. Concentration eigenvalues of K are ' ...
        'NOT eigenvalues of this orthogonal matrix.'];
    return
end
K=zeros(n);Q=zeros(n);lambda=zeros(n,1);
degree=zeros(n,1);order=zeros(n,1);kind=zeros(n,1);
% KIND=0 zonal, 1 cos(m*phi), 2 sin(m*phi).
[qnodes,qweights]=cap_quadrature(L,cap);
pos=0;
for m=0:L-1
    ell=(m:L-1)';
    h=numel(ell);
    B=zeros(L,h);
    for j=1:h
        l=ell(j);
        P=legendre(l,qnodes.');
        normfactor=sqrt((2*l+1)/2 ...
            *exp(gammaln(l-m+1)-gammaln(l+m+1)));
        B(:,j)=normfactor*P(m+1,:).';
    end
    T=B'*(qweights.*B);
    T=(T+T')/2;
    copies=1+(m>0);
    for c=1:copies
        ids=pos+(1:h);pos=pos+h;
        K(ids,ids)=T;
        degree(ids)=ell;
        order(ids)=m;
        if m>0,kind(ids)=c;end
        if cap==0||cap==pi
            V=eye(h);
            if cap==0,d=zeros(h,1);else,d=ones(h,1);end
        else
            [V,D]=eig(T);
            [d,idx]=sort(real(diag(D)),'descend');
            V=real(V(:,idx));
            for j=1:h
                [~,pivot]=max(abs(V(:,j)));
                if V(pivot,j)<0,V(:,j)=-V(:,j);end
            end
        end
        Q(ids,ids)=V;lambda(ids)=d;
    end
end
[~,perm]=sortrows([-lambda,(1:n)']);
Q=Q(:,perm);lambda=lambda(perm);
switch output
    case {'matrix','analysis'}
        U=Q.';
    case 'basis'
        U=Q;
    case 'concentration'
        U=K;
    case 'parameters'
        U=struct('concentration',K,'basis',Q,'analysis',Q.', ...
            'concentration_eigenvalues',lambda, ...
            'bandlimit',L,'cap',cap,'degree',degree, ...
            'azimuthal_order',order,'harmonic_kind',kind, ...
            'eigenvector_original_index',perm(:), ...
            'shannon_number',n*(1-cos(cap))/2, ...
            'convention','real fixed-m cos/sin, signed max-pivot positive');
end
end

function [x,w]=cap_quadrature(q,cap)
if cap==0
    x=ones(q,1);w=zeros(q,1);return
end
[d,s]=js_jacobi_recurrence(0,0,q);
J=diag(d)+diag(s,1)+diag(s,-1);
[V,D]=eig(J);
[t,ord]=sort(diag(D)/2);
v=V(1,ord).';
c=cos(cap);
x=(1+c)/2+(1-c)*t/2;
w=(1-c)*v.^2;
end
