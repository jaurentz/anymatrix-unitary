function tests=test_matrix_free
%TEST_MATRIX_FREE Forward, adjoint, multi-RHS and materialization invariants.
tests=functiontests(localfunctions);
end

function test_core_shapes_and_terminal_phases(testCase)
rngState=rng; cleanup=onCleanup(@()rng(rngState)); %#ok<NASGU>
rng(7);
for n=[1 2 3 4 7 10]
    alpha=.6*exp(1i*(1:n-1)'*.3);
    beta=exp(.47i);
    twists=mod((1:max(0,n-2))',3)==1;
    for form={'hessenberg','inverse_hessenberg','cmv','twisted'}
        rep=form{1};
        A=js_core_product(alpha,beta,rep,twists);
        O=js_core_product(alpha,beta,rep,twists,'operator');
        check_operator(testCase,A,O);
        V=opuc_matrix(n,alpha,beta,rep,twists,'operator');
        check_operator(testCase,A,V);
    end
end
end

function test_public_opuc_generator_pipeline(testCase)
n=8;beta=exp(.24i);alpha=[.2i*ones(n-1,1);beta];
mu=szego_jacobi_measure(.15,.4);
cases={
    finite_cmv(n,alpha),finite_cmv(n,alpha,'operator');
    circular_beta(n,2,12),circular_beta(n,2,12,'operator');
    jacobi_szego(n,.15,.4,'hessenberg',beta,[]),jacobi_szego(n,.15,.4,'hessenberg',beta,[],'operator');
    bernstein_szego(n,.3i,'twisted',beta,[]),bernstein_szego(n,.3i,'twisted',beta,[],'operator');
    bernstein_szego_finite(n,[.25;.2i],'cmv',beta,[]),bernstein_szego_finite(n,[.25;.2i],'cmv',beta,[],'operator');
    circular_jacobi(n,.3,.5,'cmv',beta,[]),circular_jacobi(n,.3,.5,'cmv',beta,[],'operator');
    rogers_szego(n,.5,'cmv',beta,[]),rogers_szego(n,.5,'cmv',beta,[],'operator');
    geronimus_cmv(n,.3+.2i,'cmv',beta,[]),geronimus_cmv(n,.3+.2i,'cmv',beta,[],'operator');
    periodic_cmv(n,[.2;.25i],'cmv',beta,[]),periodic_cmv(n,[.2;.25i],'cmv',beta,[],'operator');
    szego_matrix(mu,n,'inverse_hessenberg',beta,[]),szego_matrix(mu,n,'inverse_hessenberg',beta,[],'operator');
    js_measure_to_matrix(mu,n,'cmv',beta,[]),js_measure_to_matrix(mu,n,'cmv',beta,[],'operator')
    };
for j=1:size(cases,1)
    check_operator(testCase,cases{j,1},cases{j,2});
end
end

function test_fft_circulant_kicks_and_baker(testCase)
n=12;theta=.3*(0:n-1)';
check_operator(testCase,phase_circulant(n,theta),phase_circulant(n,theta,'operator'));
check_operator(testCase,cycle_propagator(n,.9),cycle_propagator(n,.9,'operator'));
check_operator(testCase,kicked_harper(n,1.1,.8),kicked_harper(n,1.1,.8,'operator'));
check_operator(testCase,kicked_rotor(n,2,.9,.3),kicked_rotor(n,2,.9,.3,'operator'));
for variant={'balazs_voros','saraceno'}
    c=variant{1};
    check_operator(testCase,quantum_baker(n,c),quantum_baker(n,c,'operator'));
end
% Independent formula for Saraceno FFT phases, including complex adjoint.
shift=.5;q=(0:n-1)'+shift;
F=exp(-2*pi*1i*(q*q')/n)/sqrt(n);
qhalf=(0:n/2-1)'+shift;
Gh=exp(-2*pi*1i*(qhalf*qhalf')/(n/2))/sqrt(n/2);
A=F'*blkdiag(Gh,Gh);
O=quantum_baker(n,'saraceno','operator');
verifyLessThan(testCase,norm(O(eye(n),'notransp')-A,'fro'),2e-11);
end

function test_floquet_circuits_and_spin_layers(testCase)
for n=[2 4]
    check_operator(testCase,dual_unitary_circuit(n,2,.34,[.2 -.3 .4 -.5],'open'), ...
        dual_unitary_circuit(n,2,.34,[.2 -.3 .4 -.5],'open','operator'));
    check_operator(testCase,xxz_floquet(n,.4,.8,'open'), ...
        xxz_floquet(n,.4,.8,'open','operator'));
end
for boundary={'open','periodic'}
    bd=boundary{1};
    check_operator(testCase,floquet_ising(4,.63,.49,bd), ...
        floquet_ising(4,.63,.49,bd,'operator'));
end
check_operator(testCase,selfdual_kicked_ising(4,[.2 -.1 .3 .45]), ...
    selfdual_kicked_ising(4,[.2 -.1 .3 .45],'operator'));
% A wrapped bond must respect ordered two-site qubit basis.
G=xxz_r_gate(.3,.8);
Q=qm_two_site_gate(4,G,4,1);
X=randn(16,2)+1i*randn(16,2);
verifyLessThan(testCase,norm(um_gate_apply(X,G,[4 1],4)-Q*X,'fro'),2e-12);
end

function test_reflections_wavelets_and_cayley(testCase)
check_operator(testCase,grover_operator(23,3),grover_operator(23,3,'operator'));
for levels=[0 1 3 5]
    check_operator(testCase,haar_wavelet(32,levels),haar_wavelet(32,levels,'operator'));
end
check_operator(testCase,cayley_jacobi(13,.2,.4,.7), ...
    cayley_jacobi(13,.2,.4,.7,'operator'));
check_operator(testCase,cayley_jacobi(1,0,0,.7), ...
    cayley_jacobi(1,0,0,.7,'operator'));
end

function test_monomial_permutations_and_diagonals(testCase)
for n=[1 2 8 19]
    check_operator(testCase,cycle_shift(n),cycle_shift(n,'operator'));
    theta=.15*(1:n)';
    check_operator(testCase,weighted_shift(n,theta),weighted_shift(n,theta,'operator'));
    check_operator(testCase,diagonal_phase(n,theta),diagonal_phase(n,theta,'operator'));
end
end

function test_larger_implicit_actions(testCase)
n=1024;
alpha=.2*exp(1i*(1:n-1)'*.1);
U=opuc_matrix(n,alpha,exp(1i*.4),'twisted',[],'operator');
x=randn(n,2)+1i*randn(n,2);
verifyEqual(testCase,U([],'size'),[n n]);
verifyLessThan(testCase,norm(U(U(x,'notransp'),'transp')-x,'fro'),1e-10);
F=phase_circulant(n,linspace(-2,2,n)','operator');
verifyLessThan(testCase,norm(F(F(x,'notransp'),'transp')-x,'fro'),1e-10);
% These orders would be impractical to materialize in circuit test loops.
C=xxz_floquet(12,.4,.8,'open','operator');
verifyEqual(testCase,C([],'size'),[4096 4096]);
v=randn(4096,1)+1i*randn(4096,1);
verifyLessThan(testCase,norm(C(C(v,'notransp'),'transp')-v),1e-10);
end

function test_bad_actions(testCase)
O=opuc_matrix(3,[.2;.3i],1,'cmv',[],'operator');
verifyError(testCase,@()O(ones(4,1),'notransp'),'um_core_apply:Dimension');
verifyError(testCase,@()O(eye(3),'bogus'),'um_core_apply:Flag');
verifyError(testCase,@()um_is_operator('bogus'),'um_is_operator:Mode');
end

function check_operator(testCase,A,O)
n=size(A,1);
verifyTrue(testCase,isa(O,'function_handle'));
verifyEqual(testCase,O([],'size'),[n n]);
X=randn(n,3)+1i*randn(n,3);
Y=randn(n,3)+1i*randn(n,3);
tol=4e-11*max(1,n);
verifyLessThan(testCase,norm(O(X,'notransp')-A*X,'fro'),tol);
verifyLessThan(testCase,norm(O(Y,'transp')-A'*Y,'fro'),tol);
verifyLessThan(testCase,norm(O(Y,'ctransp')-A'*Y,'fro'),tol);
verifyLessThan(testCase,norm(O(eye(n),'notransp')-full(A),'fro'),tol);
verifyLessThan(testCase,norm(O(O(X,'notransp'),'transp')-X,'fro'),tol);
verifyLessThan(testCase,abs(sum(conj(O(X,'notransp')).*Y,'all')- ...
    sum(conj(X).*O(Y,'transp'),'all')),tol);
end
