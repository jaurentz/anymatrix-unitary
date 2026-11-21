function [U,properties]=spherical_lune_slepian(L,phiStart,width,output)
%SPHERICAL_LUNE_SLEPIAN Exact non-axisymmetric spherical Slepian basis.
%   U=SPHERICAL_LUNE_SLEPIAN(L,PHISTART,WIDTH) returns the real
%   orthogonal ANALYSIS matrix Q' of the complete finite spherical
%   Slepian eigenbasis concentrated in the spherical LUNE
%       R={0<=theta<=pi, PHISTART<=phi<=PHISTART+WIDTH}.
%   0<=WIDTH<=2*pi; angles in radians; PHISTART default .23;
%   WIDTH default pi/3. L is the harmonic BANDLIMIT with degrees
%   ell=0,...,L-1 and N=L^2 harmonic coefficients (not matrix order).
%   N<=144 for the explicit dense reference prototype.
%
%   The real harmonics use precisely the ordering of SPHERICAL_SLEPIAN:
%     m=0: ell=0,...,L-1;
%     m=1,...,L-1: cos(m*phi), then sin(m*phi), ell=m,...,L-1.
%   Radial normalization is sqrt((2ell+1)/2*(ell-m)!/(ell+m)!)
%   times MATLAB's Condon--Shortley P_ell^m(cos(theta)).
%   Azimuthal normalization is 1/sqrt(2*pi) for m=0,
%   and 1/sqrt(pi) times cos/sin(m*phi) when m>0.
%
%   K_(a,b)=integral_R Y_a Y_b dOmega is EXACT in real arithmetic
%   from separated factors: analytic integrals of trig products on
%   [PHISTART,PHISTART+WIDTH], and L-node Gaussian quadrature in
%   x=cos(theta). For m_a+m_b even, P_la^ma*P_lb^mb is a
%   polynomial of degree <=2L-2: Gauss--Legendre is exact. For
%   m_a+m_b odd, it equals sqrt(1-x^2) times a polynomial of
%   degree <=2L-3: Gauss--Jacobi(1/2,1/2) is exact. K is a real
%   symmetric positive contraction in exact arithmetic; the
%   concentration eigenvalues ARE NOT eigenphases of U.
%
%   OUTPUT 'matrix'/'analysis': Q' (default, orthogonal)
%          'basis': Q (synthesis, orthogonal)
%          'concentration': K (NOT orthogonal/unitary)
%          'parameters': struct including K,Q,lambda,region and labels
%          'spectral': UNKNOWN eigenphases of Q', NOT lambda.
%   Empty/full lunes return the canonical Q=I, avoiding arbitrary
%   bases of repeated eigenvalues. Interior repeated eigenspaces
%   are canonically resolved by projecting coordinate vectors, with
%   a documented numerical eigenvalue tie tolerance; the sign of
%   each basis vector has its largest absolute component positive.
%
%   Setup is dense O(L^5+N^3), N=L^2; storage/action O(N^2).
%   This does NOT implement arbitrary polygons/masks, the fast
%   commuting tridiagonal cap algorithm, or any matrix-free action.
%   In contrast with approximate mask quadrature, the lune
%   integrals are exact in arithmetic; floating errors remain.
%
%   [156] Simons, Dahlen, Wieczorek, SIAM Review 48 (2006), 504--536.
%   https://doi.org/10.1137/S0036144504445765
properties={'square','unitary','orthogonal','real','normal', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','real','integer','positive','finite','<=',12});
if nargin<2||isempty(phiStart),phiStart=.23;end
if nargin<3||isempty(width),width=pi/3;end
validateattributes(phiStart,{'numeric'}, ...
    {'scalar','real','finite','>=',-1e6,'<=',1e6});
validateattributes(width,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',2*pi});
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('spherical_lune_slepian:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','analysis','basis', ...
        'concentration','parameters','spectral'}))
    error('spherical_lune_slepian:Mode','Unsupported output mode.');
end
n=L^2;
phiStart=mod(phiStart,2*pi);
if strcmp(output,'spectral')
    U=unitary_spectral_reference('spherical_lune_slepian',n);
    U.bandlimit=L;U.longitude_start=phiStart;U.longitude_width=width;
    U.reference='Simons--Dahlen--Wieczorek, SIAM Review 48 (2006)';
    U.notes=['Real orthogonal analysis of exact meridional-lune ', ...
        'concentration. Its finite matrix eigenphases are unknown. ', ...
        'Concentration eigenvalues lie in [0,1] and belong to ', ...
        'K, not to the orthogonal analysis matrix.'];
    return
end

degree=zeros(n,1);order=zeros(n,1);kind=zeros(n,1);
h=0;
for m=0:L-1
    for c=1+(m>0):1+(m>0)*2
        % c=1 for m=0 (zonal); c=2,3 for m>0 (cos,sin).
        for ell=m:L-1
            h=h+1;degree(h)=ell;order(h)=m;kind(h)=c-1;
        end
    end
end
% For m>0, the loop above runs 2:3; for m=0 it runs 1:1.
if h~=n,error('spherical_lune_slepian:Internal','Basis dimension mismatch.');end

if width==0 || width==2*pi
    if width==0,K=zeros(n);lambda=zeros(n,1);
    else,K=eye(n);lambda=ones(n,1);
    end
    Q=eye(n);perm=(1:n)';tieTolerance=0;
else
    [x0,w0]=lune_jacobi_nodes(L,0);
    [x1,w1]=lune_jacobi_nodes(L,.5);
    B0=zeros(L,n);B1=zeros(L,n);
    for ell=0:L-1
        P0=legendre(ell,x0.');
        P1=legendre(ell,x1.');
        ids=find(degree==ell);
        for k=reshape(ids,1,[])
            m=order(k);
            f=sqrt((2*ell+1)/2*exp( ...
                gammaln(ell-m+1)-gammaln(ell+m+1)));
            B0(:,k)=f*P0(m+1,:).';
            B1(:,k)=f*P1(m+1,:).';
        end
    end
    K=zeros(n);
    % Exact longitude integrals C(k)=int cos(k*phi), S(k)=int sin(k*phi).
    % Cache frequencies k in [-2L,2L], including negative indices.
    freq=(-2*L:2*L)';
    C=zeros(numel(freq),1);S=zeros(numel(freq),1);
    stop=phiStart+width;
    for k=1:numel(freq)
        q=freq(k);
        if q==0,C(k)=width;
        else
            C(k)=(sin(q*stop)-sin(q*phiStart))/q;
            S(k)=(cos(q*phiStart)-cos(q*stop))/q;
        end
    end
    for i=1:n
        m=order(i);
        if m==0,ni=1/sqrt(2*pi);
        else,ni=1/sqrt(pi);
        end
        for j=i:n
            r=order(j);
            if r==0,nj=1/sqrt(2*pi);
            else,nj=1/sqrt(pi);
            end
            if mod(m+r,2)==0
                radial=sum(w0.*B0(:,i).*B0(:,j));
            else
                radial=sum((w1./sqrt(1-x1.^2)).*B1(:,i).*B1(:,j));
            end
            if kind(i)~=2 && kind(j)~=2
                angular=(C(m-r+2*L+1)+C(m+r+2*L+1))/2;
            elseif kind(i)==2 && kind(j)==2
                angular=(C(m-r+2*L+1)-C(m+r+2*L+1))/2;
            elseif kind(i)~=2
                angular=(S(r+m+2*L+1)+S(r-m+2*L+1))/2;
            else
                angular=(S(m+r+2*L+1)+S(m-r+2*L+1))/2;
            end
            K(i,j)=ni*nj*radial*angular;
            K(j,i)=K(i,j); % exact mathematical symmetry
        end
    end
    [V,D]=eig(K);
    [lambda,perm]=sort(real(diag(D)),'descend');
    Q=real(V(:,perm));
    % Project canonical coordinate vectors onto numerically degenerate
    % eigenspaces; this avoids arbitrary eigensolver rotations.
    tieTolerance=2e-12*max(1,norm(K,2));
    i=1;
    while i<=n
        j=i;
        while j<n && abs(lambda(j)-lambda(j+1))<=tieTolerance
            j=j+1;
        end
        if j>i
            P=Q(:,i:j)*Q(:,i:j)';
            W=zeros(n,j-i+1);count=0;
            for pivot=1:n
                v=P(:,pivot);
                if count>0
                    v=v-W(:,1:count)*(W(:,1:count)'*v);
                    v=v-W(:,1:count)*(W(:,1:count)'*v);
                end
                nv=norm(v);
                if nv>1e-8
                    count=count+1;W(:,count)=v/nv;
                    if count==j-i+1,break;end
                end
            end
            if count~=j-i+1
                error('spherical_lune_slepian:Degeneracy', ...
                    'Numerical eigenspace canonicalization failed.');
            end
            Q(:,i:j)=W;
        end
        i=j+1;
    end
    for k=1:n
        [~,pivot]=max(abs(Q(:,k)));
        if Q(pivot,k)<0,Q(:,k)=-Q(:,k);end
    end
end

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
           'bandlimit',L,'longitude_start',phiStart, ...
           'longitude_width',width, ...
           'degree',degree,'azimuthal_order',order, ...
           'harmonic_kind',kind, ...
           'eigenvalue_original_index',perm(:), ...
           'eigenvalue_tie_tolerance',tieTolerance, ...
           'shannon_number',n*width/(2*pi), ...
           'quadrature','L-node exact Legendre/Jacobi-1/2; analytic phi', ...
           'convention','real Condon--Shortley m=0,cos,sin');
end
end

function [x,w]=lune_jacobi_nodes(q,alpha)
[d,s]=js_jacobi_recurrence(alpha,alpha,q);
J=diag(d/2);
if q>1,J=J+diag(s/2,1)+diag(s/2,-1);end
[V,D]=eig(J);
[x,idx]=sort(real(diag(D)));
if alpha==0,mass=2;
else,mass=pi/2;end
w=mass*real(V(1,idx).').^2;
end
