function tests=test_spherical_slepian
tests=functiontests(localfunctions);
end

function test_polar_cap_eigenbasis_and_concentration(tc)
for L=[1 2 3 7 12]
    for theta=[.3 1.1 2.3]
        par=spherical_slepian(L,theta,'parameters');
        Q=par.basis;K=par.concentration;lambda=par.concentration_eigenvalues;
        U=spherical_slepian(L,theta);
        n=L^2;
        verifyEqual(tc,size(U),[n n]);
        verifyTrue(tc,isreal(U));
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),1e-10*n);
        verifyLessThan(tc,norm(U-Q','fro'),1e-13*n);
        verifyLessThan(tc,norm(K-K','fro'),1e-12*n);
        verifyLessThan(tc,norm(Q'*K*Q-diag(lambda),'fro'),2e-10*n);
        verifyGreaterThanOrEqual(tc,min(lambda),-4e-12);
        verifyLessThanOrEqual(tc,max(lambda),1+4e-12);
        verifyLessThan(tc,abs(trace(K)-n*(1-cos(theta))/2), ...
            2e-10*n);
        verifyLessThanOrEqual(tc,max(diff(lambda)),4e-12);
        verifyEqual(tc,numel(par.degree),n);
        verifyEqual(tc,numel(par.azimuthal_order),n);
        for j=1:n
            [~,pivot]=max(abs(Q(:,j)));
            verifyGreaterThanOrEqual(tc,Q(pivot,j),-1e-12);
        end
    end
end
end

function test_independent_degree_one_cap_formula(tc)
L=2;theta=.83;c=cos(theta);
K=spherical_slepian(L,theta,'concentration');
% Physical ordering [Y00,Y10,Y11cos,Y11sin].
f=(1-c)/2;
k01=sqrt(3)*(1-c^2)/4;
k11=(1-c^3)/2;
kmm=.5-.75*c+.25*c^3;
expected=[f k01 0 0;k01 k11 0 0;0 0 kmm 0;0 0 0 kmm];
verifyLessThan(tc,norm(K-expected,'fro'),2e-13);
% The addition theorem gives an exact area-fraction Shannon number.
verifyLessThan(tc,abs(trace(K)-4*f),2e-13);
end

function test_full_empty_caps_and_deterministic_repeats(tc)
for L=[1 2 5]
    n=L^2;
    verifyLessThan(tc,norm(spherical_slepian(L,0)-eye(n),'fro'),1e-12);
    verifyLessThan(tc,norm(spherical_slepian(L,pi)-eye(n),'fro'),1e-12);
    verifyLessThan(tc,norm(spherical_slepian(L,0,'concentration'),'fro'),1e-12);
    verifyLessThan(tc,norm(spherical_slepian(L,pi,'concentration')-eye(n),'fro'),2e-12*n);
end
A=spherical_slepian(9,1.25);
B=spherical_slepian(9,1.25);
verifyEqual(tc,A,B,'AbsTol',0);
S=spherical_slepian(6,1.25,'spectral');
verifyEqual(tc,S.family,'spherical_slepian');
verifyEqual(tc,S.dimension,36);
verifyEqual(tc,S.finite.kind,'unknown');
verifyNotEmpty(tc,S.notes);
end

function test_hemisphere_complement_and_distinct_from_dpss(tc)
L=7;par=spherical_slepian(L,pi/2,'parameters');
K=par.concentration;
% South-hemisphere harmonic reflection: P_l^m(-x)=(-1)^(l+m)P_l^m(x).
D=diag((-1).^(par.degree+par.azimuthal_order));
verifyLessThan(tc,norm(K+D*K*D-eye(L^2),'fro'),3e-10);
verifyLessThan(tc,abs(trace(K)-L^2/2),3e-11);
verifyEqual(tc,numel(par.concentration_eigenvalues),L^2);
verifyError(tc,@()spherical_slepian(5,1,'operator'), ...
    'spherical_slepian:Mode');
verifyError(tc,@()spherical_slepian(5,1,'random_region'), ...
    'spherical_slepian:Mode');
end
