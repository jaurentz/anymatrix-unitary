function tests=test_spherical_lune_slepian
%TEST_SPHERICAL_LUNE_SLEPIAN Exact non-axisymmetric spherical concentration.
tests=functiontests(localfunctions);
end

function test_exact_lune_orthogonal_basis_and_concentration(tc)
for L=[1 2 3 5 9]
    n=L^2;
    for interval=[.45 1.7; 1.13 .8; 4.1 3.2]'
        start=interval(1);width=interval(2);
        P=spherical_lune_slepian(L,start,width,'parameters');
        K=P.concentration;Q=P.basis;lam=P.concentration_eigenvalues;
        U=spherical_lune_slepian(L,start,width);
        verifyEqual(tc,size(U),[n n]);
        verifyTrue(tc,isreal(U));
        verifyLessThan(tc,norm(U-Q','fro'),2e-13*n);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),2e-10*n);
        verifyLessThan(tc,norm(K-K','fro'),3e-13*n);
        verifyLessThan(tc,norm(Q'*K*Q-diag(lam),'fro'), ...
            5e-10*n);
        verifyLessThan(tc,abs(trace(K)-n*width/(2*pi)), ...
            3e-10*n);
        verifyGreaterThanOrEqual(tc,min(lam),-4e-11);
        verifyLessThanOrEqual(tc,max(lam),1+4e-11);
        verifyLessThanOrEqual(tc,max([0;diff(lam)]),3e-12);
        for j=1:n
            [~,pivot]=max(abs(Q(:,j)));
            verifyGreaterThanOrEqual(tc,Q(pivot,j),-2e-12);
        end
        verifyEqual(tc,nnz(P.degree>=0),n);
    end
end
end

function test_independent_nonaxisymmetric_double_quadrature(tc)
% This reference evaluates fully separated real harmonics on an independent
% Gauss--Legendre colatitude x Gauss--Legendre longitude tensor rule;
% neither exact-Jacobi
% parity integration nor analytic longitude formulas are reused.
for L=[2 3 5]
    for a=[.27 1.19]
        width=1.32;
        K=spherical_lune_slepian(L,a,width,'concentration');
        R=direct_double_quadrature(L,a,width,3*L+12);
        verifyLessThan(tc,norm(K-R,'fro'),2e-11*L^2);
        if L>=3
            % A lune is non-axisymmetric: several cross-m couplings exist.
            verifyGreaterThan(tc,norm(K(1:L,L+1:end),'fro'),.005);
        end
    end
end
end

function test_complements_rotation_and_boundary(tc)
for L=[1 2 5 8]
    n=L^2;start=.34;angle=1.07;
    K=spherical_lune_slepian(L,start,angle,'concentration');
    Kc=spherical_lune_slepian(L,start+angle, ...
        2*pi-angle,'concentration');
    verifyLessThan(tc,norm(K+Kc-eye(n),'fro'),4e-10*n);
    P=spherical_lune_slepian(L,start,angle,'parameters');
    D=spherical_lune_slepian(L,start+angle,angle,'parameters');
    T=eye(n);
    for m=1:L-1
        for ell=m:L-1
            idsCos=find(P.degree==ell&P.azimuthal_order==m ...
                &P.harmonic_kind==1);
            idsSin=find(P.degree==ell&P.azimuthal_order==m ...
                &P.harmonic_kind==2);
            z=cos(m*angle);s=sin(m*angle);
            T([idsCos,idsSin],[idsCos,idsSin])=[z,-s;s,z];
        end
    end
    verifyLessThan(tc,norm(D.concentration-T*K*T','fro'), ...
        3e-10*n);
    verifyLessThan(tc,norm(spherical_lune_slepian(L,start,0) ...
        -eye(n),'fro'),1e-13*n);
    verifyLessThan(tc,norm(spherical_lune_slepian(L,start,2*pi) ...
        -eye(n),'fro'),1e-13*n);
    verifyLessThan(tc,norm(spherical_lune_slepian(L,start,0, ...
        'concentration'),'fro'),1e-13*n);
    verifyLessThan(tc,norm(spherical_lune_slepian(L,start,2*pi, ...
        'concentration')-eye(n),'fro'),1e-13*n);
end
end

function test_deterministic_repeated_eigenvalue_basis_and_spectral(tc)
L=7;a=.23;w=pi/2;
P=spherical_lune_slepian(L,a,w,'parameters');
Q=spherical_lune_slepian(L,a,w,'parameters');
verifyEqual(tc,P.basis,Q.basis,'AbsTol',0);
verifyEqual(tc,P.concentration,Q.concentration,'AbsTol',0);
S=spherical_lune_slepian(L,a,w,'spectral');
verifyEqual(tc,S.family,'spherical_lune_slepian');
verifyEqual(tc,S.dimension,L^2);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
verifyGreaterThan(tc,norm(P.concentration'*P.concentration ...
    -eye(L^2),'fro'),.1); % K is not the unitary matrix.
[Z,properties]=spherical_lune_slepian();
verifyEmpty(tc,Z);
verifyTrue(tc,any(strcmp(properties,'orthogonal')));
verifyError(tc,@()spherical_lune_slepian(13), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()spherical_lune_slepian(4,.1,-.1), ...
    'MATLAB:notGreaterEqual');
verifyError(tc,@()spherical_lune_slepian(4,.1,pi,'operator'), ...
    'spherical_lune_slepian:Mode');
end

function K=direct_double_quadrature(L,phiStart,width,q)
[theta,wtheta]=independent_legendre_quadrature(q,0,pi);
x=cos(theta);wx=wtheta.*sin(theta);
[phi,wp]=independent_legendre_quadrature(q,phiStart,phiStart+width);
n=L^2;
ell=zeros(n,1);m=zeros(n,1);kind=zeros(n,1);
r=0;
for mm=0:L-1
    copies=1+(mm>0);
    for c=1:copies
        for l=mm:L-1
            r=r+1;ell(r)=l;m(r)=mm;
            if mm>0,kind(r)=c;end
        end
    end
end
B=zeros(q,n);
for j=1:n
    l=ell(j);mm=m(j);
    P=legendre(l,x.');
    f=sqrt((2*l+1)/2*exp(gammaln(l-mm+1)-gammaln(l+mm+1)));
    B(:,j)=f*P(mm+1,:).';
end
K=zeros(n);
for a=1:q
    angular=zeros(1,n);
    for j=1:n
        if m(j)==0,angular(j)=1/sqrt(2*pi);
        elseif kind(j)==1,angular(j)=cos(m(j)*phi(a))/sqrt(pi);
        else,angular(j)=sin(m(j)*phi(a))/sqrt(pi);
        end
    end
    V=B.*angular;
    K=K+wp(a)*(V'*(wx.*V));
end
end

function [x,w]=independent_legendre_quadrature(q,a,b)
% Golub--Welsch with classical Legendre symmetric tridiagonal.
k=(1:q-1)';
off=k./sqrt(4*k.^2-1);
J=diag(off,1)+diag(off,-1);
[V,D]=eig(J);
[z,idx]=sort(diag(D));
x=(a+b)/2+(b-a)*z/2;
w=(b-a)*V(1,idx).'.^2;
end
