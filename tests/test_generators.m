function tests=test_generators
%TEST_GENERATORS Matrix generator tests and analytic reference spectra.
tests=functiontests(localfunctions);
end

function test_cycle_shift_and_phase(testCase)
for n=[1 2 3 4 8 17]
    P=cycle_shift(n);
    assertUnitary(testCase,P);
    verifyTrue(testCase,issparse(P));
    verifyLessThan(testCase,norm(full(P^n-speye(n)),'fro'),1e-12);
    for k=1:n-1
        verifyLessThan(testCase,abs(trace(P^k)),1e-12);
    end
    theta=linspace(-2,2,n).';
    D=diagonal_phase(n,theta);
    assertUnitary(testCase,D);
    verifyLessThan(testCase,norm(diag(D)-exp(1i*theta)),1e-13);
    verifyLessThan(testCase,norm(full(D-diag(diag(D))),'fro'),1e-13);
end
end

function test_cycle_propagator(testCase)
for n=[1 2 3 5 8 17]
    for t=[0 0.3 2.1]
        P=cycle_shift(n);
        U=cycle_propagator(n,t);
        assertUnitary(testCase,U);
        verifyLessThan(testCase,norm(U-expm(-1i*t*full(P+P')),'fro'),3e-12*n);
        lambda=exp(-2i*t*cos(2*pi*(0:n-1)'/n));
        verifyLessThan(testCase,abs(trace(U)-sum(lambda)),3e-12*n);
    end
end
end

function test_cmv_is_shared_implementation(testCase)
for n=[1 2 3 7 12]
    alpha=[.45*exp(1i*(1:n-1)'); exp(1i*.13)];
    C=finite_cmv(n,alpha);
    G=opuc_matrix(n,alpha(1:end-1),alpha(end),'cmv');
    verifyEqual(testCase,C,G,'AbsTol',1e-13);
    verifyTrue(testCase,issparse(C));
    assertUnitary(testCase,C);
end
end

function test_quantum_walk(testCase)
for n=[1 2 3 5 8]
    for angles=[.3 -.7; 0 0; pi/3 pi/9]'
        U=split_step_walk(n,angles(1),angles(2));
        verifyEqual(testCase,size(U),[2*n 2*n]);
        verifyTrue(testCase,isreal(U));
        verifyTrue(testCase,issparse(U));
        assertUnitary(testCase,U);
    end
end
end

function test_haar_unitary(testCase)
for n=[1 2 4 11]
    state=rng;
    U=haar_unitary(n,123);
    verifyEqual(testCase,rng,state);
    assertUnitary(testCase,U);
    verifyEqual(testCase,U,haar_unitary(n,123),'AbsTol',0);
end
end

function test_jacobi_public_api(testCase)
mu=szego_jacobi_measure(0,0);
for n=[1 2 8 12]
    [C,alpha]=szego_matrix(mu,n,'cmv',1);
    J=jacobi_szego(n,0,0,'cmv',1);
    verifyEqual(testCase,C,J,'AbsTol',1e-13);
    verifyEqual(testCase,alpha,szego_verblunsky(mu,n-1),'AbsTol',1e-13);
    assertUnitary(testCase,J);
end
end

function test_anymatrix_metadata(testCase)
fun={@cycle_shift,@diagonal_phase,@cycle_propagator,@finite_cmv, ...
     @split_step_walk,@haar_unitary,@jacobi_szego,@opuc_matrix};
for k=1:numel(fun)
    [U,prop]=fun{k}();
    verifyEmpty(testCase,U);
    verifyTrue(testCase,iscell(prop));
    verifyTrue(testCase,any(strcmp(prop,'unitary')));
    verifyTrue(testCase,any(strcmp(prop,'square')));
end
end

function assertUnitary(testCase,U)
n=size(U,1);
verifyEqual(testCase,size(U,2),n);
verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),3e-12*max(n,1));
end
