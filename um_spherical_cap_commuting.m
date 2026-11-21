function C=um_spherical_cap_commuting(L,cap)
%UM_SPHERICAL_CAP_COMMUTING Exact Grünbaum tridiagonal Slepian blocks.
%   For maximal harmonic degree L-1 and fixed order m:
%     T(l,l)=cos(cap)*l*(l+1),
%     T(l,l+1)=(L^2-(l+1)^2)
%         *sqrt(((l+1)^2-m^2)/((2*l+1)*(2*l+3))).
%   T is -1 times the Grünbaum operator matrix in Simons et al.
%   (SIAM Rev. 48, 2006, section 5.4), with their maximal degree
%   parameter set to L-1. [156]
%   https://doi.org/10.1137/S0036144504445765
%
%   T commutes with the fixed-m polar-cap concentration K. Its
%   eigenvectors are stably identifiable even when eigenvalues of K
%   cluster near zero or one; the monotone T eigenvalues must not be
%   confused with concentration factors or analysis eigenphases.
%   Finite lambda are computed by Gaussian cap quadrature as Rayleigh
%   quotients of the T eigenvectors (rather than eig(K)).
%
%   Block setup is O(L^4) using generic dense eig(T), O(L^3)
%   stored eigenvectors, and O(L^3) matvec per right-hand side.
%   No claim of a specialized tridiagonal O(L^3) TOTAL eigensolver,
%   a fast spherical harmonic transform or approximate Krylov unitarity.
n=L^2;
c=cos(cap);
[ds,es]=js_jacobi_recurrence(0,0,L);
J=diag(ds/2);
if L>1,J=J+diag(es/2,1)+diag(es/2,-1);end
[Vq,Dq]=eig(J);
[r,idx]=sort(real(diag(Dq)));
x=(1+c)/2+(1-c)*r/2;
w=(1-c)*real(Vq(1,idx).').^2;
blocks=cell(2*L-1,1);
tridiagonals=cell(L,1);
localEigen=zeros(n,1);
preLambda=zeros(n,1);
degree=zeros(n,1);order=zeros(n,1);kind=zeros(n,1);
pos=0;count=0;
for m=0:L-1
    ell=(m:L-1)';
    h=numel(ell);
    diagonal=c*ell.*(ell+1);
    j=(m:L-2)';
    off=(L^2-(j+1).^2).*sqrt( ...
        ((j+1).^2-m^2)./((2*j+1).*(2*j+3)));
    T=diag(diagonal);
    if h>1,T=T+diag(off,1)+diag(off,-1);end
    tridiagonals{m+1}=sparse(T);
    if cap==0 || cap==pi
        V=eye(h);
        concentration=double(cap==pi)*ones(h,1);
        chi=diag(T);
    else
        [V,D]=eig(T);
        [chi,idx]=sort(real(diag(D)),'descend');
        V=real(V(:,idx));
        for k=1:h
            [~,pivot]=max(abs(V(:,k)));
            if V(pivot,k)<0,V(:,k)=-V(:,k);end
        end
        B=zeros(L,h);
        for k=1:h
            l=ell(k);
            P=legendre(l,x.');
            f=sqrt((2*l+1)/2*exp( ...
                gammaln(l-m+1)-gammaln(l+m+1)));
            B(:,k)=f*P(m+1,:).';
        end
        BV=B*V;
        concentration=sum(w.*(BV.^2),1).';
    end
    copies=1+(m>0);
    for kindIndex=1:copies
        count=count+1;
        ids=pos+(1:h);
        pos=pos+h;
        blocks{count}=struct('ids',ids,'V',V,'order',m, ...
            'kind',kindIndex*(m>0), ...
            'grunbaum_eigenvalues',chi);
        preLambda(ids)=concentration;
        localEigen(ids)=chi;
        degree(ids)=ell;
        order(ids)=m;
        if m>0,kind(ids)=kindIndex;end
    end
end
if pos~=n,error('um_spherical_cap_commuting:Internal','Unexpected basis order.');end
[~,perm]=sortrows([-preLambda,(1:n)']);
C=struct('blocks',{blocks},'tridiagonal_blocks',{tridiagonals}, ...
    'concentration_eigenvalues',preLambda(perm), ...
    'grunbaum_eigenvalues',localEigen(perm), ...
    'sorted_permutation',perm(:), ...
    'degree',degree,'azimuthal_order',order, ...
    'harmonic_kind',kind,'bandlimit',L,'cap',cap, ...
    'dimension',n,'shannon_number',n*(1-c)/2, ...
    'notes',['Fixed-order exact commuting tridiagonal; its ', ...
    'eigenvalues are NOT concentration factors. Local Rayleigh ', ...
    'factors use cap quadrature; global basis sorted by factors. ', ...
    'Dense eigenproblem per tridiagonal block, no big K formed.']);
end
