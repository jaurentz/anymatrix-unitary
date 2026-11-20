function tests=test_geometric_fem_wave
%TEST_GEOMETRIC_FEM_WAVE Dirichlet finite-element wave propagator verification.
tests=functiontests(localfunctions);
end

function test_exact_finite_element_frequencies_and_metric(tc)
for n=[1,2,6,11]
    dx=2/(n+1);h=.08*dx;
    for method={'midpoint','gauss2','gauss3','verlet'}
        name=method{1};
        U=geometric_fem_wave(n,2,1.4,.6,h,name);
        P=geometric_fem_wave(n,2,1.4,.6,h,name,'pencil');
        spec=geometric_fem_wave(n,2,1.4,.6,h,name,'spectral');
        specEig=spec.finite.eigenvalues();
        verifyEqual(tc,size(U),[2*n,2*n]);
        verifyEqual(tc,spec.dimension,2*n);
        verifyEqual(tc,spec.finite.kind,'exact');
        verifyLessThan(tc,abs(spec.finite.moments(0)-1),1e-13);
        verifyLessThan(tc,norm(U'*U-eye(2*n),'fro'),5e-10*n);
        J=[zeros(n),eye(n);-eye(n),zeros(n)];
        verifyLessThan(tc,norm(U'*J*U-J,'fro'),5e-10*n);
        verifyLessThan(tc,norm(P.M'*J*P.M-J,'fro'),1e-8*n);
        verifyGreaterThan(tc,min(eig(P.G)),0);
        verifyLessThan(tc,norm(P.M'*P.G*P.M-P.G,'fro')/ ...
            norm(P.G,'fro'),1e-8);
        verifyLessThan(tc,norm(P.A'* (P.B\P.A)-P.B,'fro')/ ...
            norm(P.B,'fro'),1e-8);
        verifyLessThan(tc,norm(P.T*P.M/P.T-U,'fro'),2e-9*n);
        verifyLessThan(tc,norm(P.T'*P.T-P.G,'fro')/ ...
            norm(P.G,'fro'),1e-10);
        w2=sort(eig(full(P.stiffness),full(P.mass)));
        verifyLessThan(tc,norm(sort(spec.frequencies.^2)-w2), ...
            2e-10*max(1,max(w2)));
        eigs=eig(U);
        for j=1:2*n
            verifyLessThan(tc,min(abs(eigs-specEig(j))),2e-8);
        end
    end
end
end

function test_verlet_finite_element_palindromic_qep(tc)
for n=[1,2,8]
    dx=1/(n+1);h=.1*dx;
    R=geometric_fem_wave(n,1,1,.5,h,'verlet','polynomial');
    verifyLessThan(tc,norm(full(R.A0-R.A2),'fro'),1e-14);
    verifyLessThan(tc,norm(full(R.A1-R.A1'),'fro'),1e-14);
    verifyLessThan(tc,norm(full(R.A0-R.A0'),'fro'),1e-14);
    C=[zeros(n),eye(n);-eye(n),-full(R.A0\R.A1)];
    roots=eig(C);
    S=geometric_fem_wave(n,1,1,.5,h,'verlet','spectral');
    SEig=S.finite.eigenvalues();
    for j=1:2*n
        verifyLessThan(tc,min(abs(roots-SEig(j))),1e-7);
    end
    P=geometric_fem_wave(n,1,1,.5,h,'verlet','pencil');
    Mass=full(P.mass);K=full(P.stiffness);
    velocityStep=[eye(n)-.5*h*h*(Mass\K),h*(Mass\eye(n)); ...
        -h*K+.25*h^3*K*(Mass\K),eye(n)-.5*h*h*(K/Mass)];
    verifyLessThan(tc,norm(P.M-velocityStep,'fro'),2e-9*n);
end
end

function test_fast_dst_operator(tc)
for n=[1 2 7 11]
    h=.1/(n+1);
    for method={'midpoint','gauss2','gauss3','verlet'}
        name=method{1};
        U=geometric_fem_wave(n,1,1,.8,h,name);
        O=geometric_fem_wave(n,1,1,.8,h,name,'operator');
        verifyTrue(tc,isa(O,'function_handle'));
        verifyEqual(tc,O([],'size'),[2*n,2*n]);
        X=randn(2*n,3)+1i*randn(2*n,3);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'), ...
            5e-11*n);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            5e-11*n);
        verifyLessThan(tc,norm(O(X,'transp')-U'*X,'fro'), ...
            5e-11*n);
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
            5e-11*n);
        verifyLessThan(tc,norm(O(speye(2*n),'notransp')-U,'fro'), ...
            5e-11*n);
    end
end
end

function test_large_fast_dst_operator(tc)
n=2048;
O=geometric_fem_wave(n,1,1,.8,[], ...
    'gauss3','operator');
X=randn(2*n,2)+1i*randn(2*n,2);
verifyEqual(tc,O([],'size'),[2*n,2*n]);
verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
    3e-10);
verifyError(tc,@()O(zeros(2*n-1,1),'notransp'), ...
    'um_geometric_wave_apply:Dimension');
verifyError(tc,@()O(X,'foo'), ...
    'um_geometric_wave_apply:Flag');
end

function test_validations(tc)
[u,p]=geometric_fem_wave();
verifyEmpty(tc,u);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyEqual(tc,size(geometric_fem_wave(4)),[8 8]);
verifyError(tc,@()geometric_fem_wave(5,1,1,0,1,'verlet'), ...
    'geometric_fem_wave:Stability');
verifyError(tc,@()geometric_fem_wave(5,1,1,0,.1,'gauss2','polynomial'), ...
    'geometric_fem_wave:Polynomial');
verifyError(tc,@()geometric_fem_wave(5,1,1,0,.1,'foo'), ...
    'geometric_fem_wave:Method');
end
