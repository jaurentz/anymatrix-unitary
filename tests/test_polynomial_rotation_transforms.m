function tests=test_polynomial_rotation_transforms
%TEST_POLYNOMIAL_ROTATION_TRANSFORMS Mathematical identities for new families.
tests=functiontests(localfunctions);
end

function test_hermite_quadrature_transform(tc)
for n=[1 2 3 8 21]
    Q=gauss_hermite_transform(n);
    orthogonal(tc,Q,6e-12*n);
    J=diag(sqrt((1:n-1)/2),1)+diag(sqrt((1:n-1)/2),-1);
    if n==1,J=0;end
    nodes=sort(eig(J));
    verifyLessThan(tc,norm(Q*J*Q'-diag(nodes),'fro'),8e-12*n);
    verifyTrue(tc,all(Q(:,1)>=-1e-15));
    if n>1
        verifyLessThan(tc,norm(Q(:,2)-sqrt(2)*nodes.*Q(:,1)), ...
            1e-11*n);
    end
end
end

function test_laguerre_quadrature_transform(tc)
for n=[1 2 6 15]
    for a=[0 .5 2.5]
        Q=gauss_laguerre_transform(n,a);
        orthogonal(tc,Q,5e-12*n);
        d=2*(0:n-1)+a+1;
        e=sqrt((1:n-1).*((1:n-1)+a));
        J=diag(d)+diag(e,1)+diag(e,-1);
        nodes=sort(eig(J));
        verifyGreaterThan(tc,min(nodes),0);
        verifyLessThan(tc,norm(Q*J*Q'-diag(nodes),'fro'),5e-11*n);
        if n>1
            verifyLessThan(tc,norm(Q(:,2) - ...
                ((nodes-(a+1))/sqrt(a+1)).*Q(:,1)),1e-11*n);
        end
    end
end
verifyEqual(tc,gauss_laguerre_transform(6), ...
    gauss_laguerre_transform(6,0),'AbsTol',1e-12);
end

function test_associated_legendre_fixed_order(tc)
for n=[1 2 5 14]
    for m=[0 1 3 8]
        Q=associated_legendre_transform(n,m);
        orthogonal(tc,Q,5e-12*n);
        k=1:n-1;
        e=sqrt(k.*(k+2*m)./((2*k+2*m).^2-1));
        J=diag(e,1)+diag(e,-1);
        if n==1,J=0;end
        nodes=sort(eig(J));
        verifyLessThan(tc,norm(Q*J*Q'-diag(nodes),'fro'),7e-12*n);
        verifyLessThan(tc,max(abs(nodes+flipud(nodes))),3e-12*n);
        if n>1
            verifyLessThan(tc,norm(Q(:,2)-(nodes/e(1)).*Q(:,1)), ...
                2e-11*n);
        end
    end
end
verifyEqual(tc,associated_legendre_transform(8,0), ...
    associated_legendre_transform(8),'AbsTol',1e-12);
end

function test_wigner_representation_and_exact_character(tc)
angles=[.3 .7 .4; .15 .66 -.48; 0 .2 0];
for ell=[0 1 2 4 8]
    m=(-ell:ell)';
    for j=1:size(angles,1)
        a=angles(j,1);b=angles(j,2);g=angles(j,3);
        D=wigner_rotation(ell,a,b,g);
        unitary(tc,D,2e-11*max(1,ell));
        delta=2*acos(cos(b/2)*cos((a+g)/2));
        for k=1:4
            verifyLessThan(tc,abs(trace(D^k)-sum(exp(-1i*k*m*delta))), ...
                4e-10*max(1,ell));
        end
        verifyLessThan(tc,abs(det(D)-1),3e-10*max(1,ell));
    end
end
end

function test_spherical_rotation_blocks_and_moments(tc)
a=.3;b=.7;g=.4;
delta=2*acos(cos(b/2)*cos((a+g)/2));
for L=[1 2 3 5 8]
    U=spherical_harmonic_rotation(L,a,b,g);
    verifyTrue(tc,issparse(U));
    verifyEqual(tc,size(U),[L^2 L^2]);
    unitary(tc,U,8e-11*L^2);
    offset=0;
    for ell=0:L-1
        idx=offset+(1:2*ell+1);
        verifyLessThan(tc,norm(full(U(idx,idx))- ...
            wigner_rotation(ell,a,b,g),'fro'),1e-12);
        offset=offset+2*ell+1;
    end
    for k=1:4
        expected=abs(sum(exp(-1i*k*(0:L-1)'*delta)))^2;
        verifyLessThan(tc,abs(trace(U^k)-expected),5e-10*L^2);
    end
end
end

function test_periodic_sinc_shift_and_fast_action(tc)
for n=[1 3 5 9 17]
    for delta=[0 .375 1 -1.25]
        Q=periodic_sinc_shift(n,delta);
        orthogonal(tc,Q,5e-11*n);
        verifyTrue(tc,isreal(Q));
        F=periodic_sinc_shift(n,delta,'operator');
        verifyEqual(tc,F([],'size'),[n n]);
        X=randn(n,2)+1i*randn(n,2);
        verifyLessThan(tc,norm(F(X,'notransp')-Q*X,'fro'), ...
            4e-11*n);
        verifyLessThan(tc,norm(F(X,'transp')-Q'*X,'fro'), ...
            4e-11*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-Q'*X,'fro'), ...
            4e-11*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'), ...
            4e-11*n);
    end
    A=periodic_sinc_shift(n,.3);
    B=periodic_sinc_shift(n,.4);
    C=periodic_sinc_shift(n,.7);
    verifyLessThan(tc,norm(A*B-C,'fro'),5e-11*n);
    P=circshift(eye(n),1,1);
    verifyLessThan(tc,norm(periodic_sinc_shift(n,1)-P,'fro'), ...
        5e-11*n);
end
verifyError(tc,@()periodic_sinc_shift(4,.3), ...
    'periodic_sinc_shift:Order');
F=periodic_sinc_shift(5,.2,'operator');
verifyError(tc,@()F(zeros(4,2),'notransp'), ...
    'um_periodic_sinc_apply:Dimension');
verifyError(tc,@()F(zeros(5,2),'bad'), ...
    'um_periodic_sinc_apply:Flag');
end

function orthogonal(tc,Q,tol)
verifyTrue(tc,isreal(Q));
unitary(tc,Q,tol);
end

function unitary(tc,U,tol)
n=size(U,1);
verifyEqual(tc,size(U),[n n]);
verifyLessThan(tc,norm(full(U'*U-speye(n)),'fro'),tol);
end

function test_exact_spectral_reference_interfaces(tc)
a=.3;b=.7;g=.4;
for ell=[0 1 3 7]
    S=wigner_rotation(ell,a,b,g,'spectral');
    D=wigner_rotation(ell,a,b,g);
    verifyEqual(tc,S.dimension,2*ell+1);
    verifyEqual(tc,S.finite.kind,'exact');
    z=S.finite.eigenvalues();
    verifyEqual(tc,numel(z),2*ell+1);
    verifyLessThan(tc,abs(sum(z.^3)-trace(D^3)),2e-10*(2*ell+1));
    verifyLessThan(tc,abs(S.finite.moments(3)-trace(D^3)/(2*ell+1)), ...
        2e-10*(2*ell+1));
end
for L=[1 2 4 7]
    S=spherical_harmonic_rotation(L,a,b,g,'spectral');
    U=spherical_harmonic_rotation(L,a,b,g);
    verifyEqual(tc,S.dimension,L^2);
    z=S.finite.eigenvalues();
    verifyEqual(tc,numel(z),L^2);
    for k=0:4
        verifyLessThan(tc,abs(sum(z.^k)-trace(U^k)),3e-10*L^2);
        verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/L^2), ...
            3e-10*L^2);
    end
end
for n=[1 3 5 9]
    S=periodic_sinc_shift(n,.37,'spectral');
    U=periodic_sinc_shift(n,.37);
    verifyEqual(tc,S.dimension,n);
    z=S.finite.eigenvalues();
    for k=1:3
        verifyLessThan(tc,abs(sum(z.^k)-trace(U^k)),3e-10*n);
    end
end
end
