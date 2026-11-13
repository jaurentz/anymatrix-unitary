function tests=test_structural_spectral_identities
%TEST_STRUCTURAL_SPECTRAL_IDENTITIES Algebraic and exact-spectrum checks for matrix families.
tests=functiontests(localfunctions);
end

function test_phase_circulant(testCase)
for n=[1 2 3 8 15]
    theta=(1:n)'*0.43+sin((1:n)')*0.2;
    U=phase_circulant(n,theta);
    verify_unitarity(testCase,U,1e-11);
    for k=1:min(n,5)
        verifyLessThan(testCase,abs(trace(U^k)-sum(exp(1i*k*theta))),1e-10*n);
    end
    verifyLessThan(testCase,norm(U(:,2:end)-circshift(U(:,1:end-1),1,1),'fro'),1e-11*n);
end
end

function test_weighted_shift(testCase)
for n=[1 2 3 9]
    theta=0.1+(1:n)'/7;
    U=weighted_shift(n,theta);
    verifyTrue(testCase,issparse(U));
    verifyEqual(testCase,nnz(U),n);
    verify_unitarity(testCase,U,1e-12);
    verifyLessThan(testCase,norm(full(U^n-exp(1i*sum(theta))*speye(n)),'fro'),1e-10);
end
end

function test_orthogonal_rotation_spectra(testCase)
for n=[1 2 3 8 11]
    angles=0.14*(1:floor(n/2))';
    for epsilon=[-1 1]
        U=orthogonal_blocks(n,angles,epsilon);
        verifyTrue(testCase,isreal(U));
        verifyTrue(testCase,issparse(U));
        verify_unitarity(testCase,U,1e-12);
        expected=[exp(1i*angles);exp(-1i*angles)];
        if mod(n,2),expected=[expected;epsilon];end
        for k=1:min(n,5)
            verifyLessThan(testCase,abs(trace(U^k)-sum(expected.^k)),1e-11*n);
        end
    end
end
end

function test_real_haar(testCase)
s0=rng;
for n=[1 2 5 13]
    U=haar_orthogonal(n,123);
    verifyTrue(testCase,isreal(U));
    verify_unitarity(testCase,U,3e-12);
    verifyLessThan(testCase,norm(U-haar_orthogonal(n,123),'fro'),1e-14);
end
verifyEqual(testCase,rng,s0);
end

function test_dirac_automaton_dispersion(testCase)
for n=[1 2 3 8 11]
    for m=[0 0.4 0.8 1]
        U=dirac_automaton(n,m);
        verifyEqual(testCase,size(U),[2*n,2*n]);
        verifyTrue(testCase,issparse(U));
        verify_unitarity(testCase,U,5e-12);
        nu=sqrt(1-m^2);
        theta=(0:n-1)'*2*pi/n;
        omega=acos(nu*cos(theta));
        reference=[exp(1i*omega);exp(-1i*omega)];
        for k=1:min(2*n,6)
            verifyLessThan(testCase,abs(trace(U^k)-sum(reference.^k)),3e-10*n);
        end
        if m>0 && m<1
            verifyLessThanOrEqual(testCase,max(abs(real(eig(full(U))))),nu+1e-12);
        end
    end
end
end

function test_grover_exact_multiplicities(testCase)
for n=[2 3 6 13]
    for m=1:n-1
        U=grover_operator(n,m);
        verifyTrue(testCase,isreal(U));
        verify_unitarity(testCase,U,1e-11);
        angle=2*asin(sqrt(m/n));
        reference=[exp(1i*angle);exp(-1i*angle);ones(m-1,1);-ones(n-m-1,1)];
        for k=1:min(n,6)
            verifyLessThan(testCase,abs(trace(U^k)-sum(reference.^k)),3e-10*n);
        end
    end
end
end

function test_quantum_baker(testCase)
for n=[2 4 6 10 20]
    for variant={'balazs_voros','saraceno'}
        U=quantum_baker(n,variant{1});
        verify_unitarity(testCase,U,2e-11);
    end
end
verifyError(testCase,@()quantum_baker(3),'quantum_baker:EvenOrder');
verifyError(testCase,@()quantum_baker(4,'other'),'quantum_baker:Variant');
end

function test_kicked_harper_degenerate_cases(testCase)
for n=[1 2 7 16]
    U=kicked_harper(n,0,0);
    verifyLessThan(testCase,norm(U-eye(n),'fro'),1e-11);
    verify_unitarity(testCase,kicked_harper(n,.4,2.1),2e-11);
    theta=2*pi*(0:n-1)'/n;
    L=0.71;
    U=kicked_harper(n,0,L);
    verifyLessThan(testCase,norm(U-diag(exp(-1i*L*cos(theta))),'fro'),2e-11);
end
end

function test_julia_operator(testCase)
for n=[1 2 3 8]
    T=0.35*full(cycle_shift(n));
    U=julia_unitary(n,T);
    verifyEqual(testCase,size(U),[2*n,2*n]);
    verify_unitarity(testCase,U,3e-11);
    verifyLessThan(testCase,norm(U(1:n,1:n)-T,'fro'),1e-12);
    Z=julia_unitary(n,zeros(n));
    verifyLessThan(testCase,norm(Z-[zeros(n),eye(n);eye(n),zeros(n)],'fro'),1e-11);
    T=(0.1+0.1i)*full(cycle_shift(n));
    verify_unitarity(testCase,julia_unitary(n,T),3e-11);
end
verifyError(testCase,@()julia_unitary(3,1.1*eye(3)),'julia_unitary:Contraction');
end

function test_cayley_of_jacobi(testCase)
for n=[1 2 4 12]
    for ab=[0 0;0.5 -0.3;2 1]'
        a=ab(1); b=ab(2);
        U=cayley_jacobi(n,a,b,.8);
        verify_unitarity(testCase,U,5e-11);
        [d,s]=js_jacobi_recurrence(a,b,n);
        J=diag(d)+diag(s,1)+diag(s,-1);
        x=eig(J);
        lambda=(x-.8i)./(x+.8i);
        for k=1:min(n,6)
            verifyLessThan(testCase,abs(trace(U^k)-sum(lambda.^k)),2e-10*n);
        end
    end
end
end

function test_polar_factor(testCase)
for n=[1 2 5 9]
    A=eye(n)+.4*diag(ones(n-1,1),1);
    U=polar_unitary(n,A);
    verifyTrue(testCase,isreal(U));
    verify_unitarity(testCase,U,5e-11);
    H=U'*A;
    verifyLessThan(testCase,norm(H-H','fro'),1e-11);
    verifyGreaterThanOrEqual(testCase,min(eig((H+H')/2)),-1e-11);
    B=A+1i*.12*eye(n);
    verify_unitarity(testCase,polar_unitary(n,B),5e-11);
end
verifyError(testCase,@()polar_unitary(2,zeros(2)),'polar_unitary:Singular');
end

function test_haar_wavelet(testCase)
for n=[1 2 4 8 32]
    for levels=0:round(log2(n))
        U=haar_wavelet(n,levels);
        verifyTrue(testCase,issparse(U));
        verifyTrue(testCase,isreal(U));
        verify_unitarity(testCase,U,1e-11);
        if levels==0
            verifyLessThan(testCase,norm(full(U)-eye(n),'fro'),1e-12);
        end
    end
    U=haar_wavelet(n);
    c=full(U*ones(n,1));
    verifyLessThan(testCase,norm(c-[sqrt(n);zeros(n-1,1)]),1e-11);
end
verifyError(testCase,@()haar_wavelet(5),'haar_wavelet:Order');
end

function verify_unitarity(testCase,U,tol)
n=size(U,1);
verifyEqual(testCase,size(U),[n n]);
verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),tol);
end
