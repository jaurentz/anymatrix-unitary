function tests=test_hadamard_rank1_etna
%TEST_HADAMARD_RANK1_ETNA Independent tests for citation-led families.
tests=functiontests(localfunctions);
end

function test_complex_hadamard_modulus_and_orthogonality(testCase)
for n=[1 2 3 4 6 8 9 10 12 15 16 18]
    U=complex_hadamard(n);
    verifyEqual(testCase,size(U),[n,n]);
    verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-11*n);
    verifyLessThan(testCase,max(abs(abs(U(:))-1/sqrt(n))),2e-12*n);
    F=complex_hadamard(n,[],[],'operator');
    verifyEqual(testCase,F([],'size'),[n n]);
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),2e-11*n);
    verifyLessThan(testCase,norm(F(X,'transp')-U'*X,'fro'),2e-11*n);
    verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),2e-11*n);
end
end

function test_hadamard_dita_formula_and_fourier_limit(testCase)
n=12;m=3;r=n/m;
twists=(1:r)'*(1:m)*.117;
U=complex_hadamard(n,m,twists);
Fm=discrete_fourier(m);
Fr=discrete_fourier(r);
T=zeros(n,n);
for j=1:m
    idx=(j-1)*r+(1:r);
    T(idx,idx)=diag(exp(1i*twists(:,j)))*Fr;
end
expected=kron(Fm,eye(r))*T;
verifyLessThan(testCase,norm(U-expected,'fro'),2e-12*n);
verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-12*n);
for n=[1 2 3 7 11]
    U=complex_hadamard(n,1,zeros(n,1));
    verifyLessThan(testCase,norm(U-discrete_fourier(n),'fro'),2e-12*n);
    S=complex_hadamard(n,1,zeros(n,1),'spectral');
    verifyEqual(testCase,S.finite.kind,'exact');
    verifyLessThan(testCase,abs(sum(S.finite.eigenvalues())-trace(U)),2e-10*n);
end
S=complex_hadamard(12,3,twists,'spectral');
verifyEqual(testCase,S.finite.kind,'unknown');
end

function test_hadamard_large_matrix_free(testCase)
n=8192;
F=complex_hadamard(n,128,[],'operator');
X=randn(n,2)+1i*randn(n,2);
Y=F(X,'notransp');
verifyLessThan(testCase,norm(F(Y,'transp')-X,'fro'),1e-9);
end

function test_hadamard_input_validation(testCase)
verifyError(testCase,@()complex_hadamard(10,3),'complex_hadamard:Factor');
verifyError(testCase,@()complex_hadamard(6,2,zeros(2,3)),'complex_hadamard:Twists');
verifyError(testCase,@()complex_hadamard(6,2,1i*ones(3,2)),'complex_hadamard:Twists');
verifyError(testCase,@()complex_hadamard(6,2,[],'invalid'),'complex_hadamard:Mode');
F=complex_hadamard(6,2,[],'operator');
verifyError(testCase,@()F(zeros(5,1),'notransp'),'um_hadamard_apply:Dimension');
verifyError(testCase,@()F(zeros(6,1),'invalid'),'um_hadamard_apply:Flag');
end

function test_rank_one_hessenberg_and_determinant(testCase)
for n=[1 2 3 8 14]
    alpha=.25*exp(1i*(1:n-1)'*.127);
    beta=exp(.31i);
    H=opuc_matrix(n,alpha,beta,'hessenberg');
    for phi=[-pi -.4 0 1.1]
        U=unitary_rank_one_update(n,phi,alpha,beta);
        tau=exp(1i*phi)-1;
        R=eye(n);R(n,n)=exp(1i*phi);
        verifyLessThan(testCase,norm(U-H*R,'fro'),2e-11*n);
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-11*n);
        verifyLessThan(testCase,norm(tril(U,-2),'fro'),2e-11*n);
        verifyLessThan(testCase,abs(det(U)/det(H)-exp(1i*phi)),2e-9*n);
        verifyLessThan(testCase,norm(U-H-tau*(H(:,n))*[zeros(1,n-1) 1],'fro'), ...
            3e-11*n);
    end
end
end

function test_rank_one_general_vector_and_implicit(testCase)
n=11;alpha=.1*exp(.2i*(1:n-1)');beta=exp(.33i);
v=(1:n)'+1i*(n:-1:1)';
phi=.91;
U=unitary_rank_one_update(n,phi,alpha,beta,v);
H=opuc_matrix(n,alpha,beta,'hessenberg');
v=v/norm(v);
verifyLessThan(testCase,norm(U-H-(exp(1i*phi)-1)*(H*v)*v','fro'),2e-11*n);
verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-11*n);
F=unitary_rank_one_update(n,phi,alpha,beta,v,'operator');
X=randn(n,3)+1i*randn(n,3);
verifyEqual(testCase,F([],'size'),[n n]);
verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),2e-11*n);
verifyLessThan(testCase,norm(F(X,'transp')-U'*X,'fro'),2e-11*n);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-11*n);
S=unitary_rank_one_update(n,phi,alpha,beta,v,'spectral');
verifyEqual(testCase,S.finite.kind,'unknown');
end

function test_rank_one_cyclic_exact_spectrum(testCase)
for n=[1 2 3 7 14]
    beta=exp(.37i);phi=.88;
    U=unitary_rank_one_update(n,phi,zeros(n-1,1),beta);
    S=unitary_rank_one_update(n,phi,zeros(n-1,1),beta,[],'spectral');
    z=S.finite.eigenvalues();
    verifyEqual(testCase,S.finite.kind,'exact');
    verifyEqual(testCase,numel(z),n);
    verifyLessThan(testCase,max(arrayfun(@(x)min(abs(x-z)),eig(U))),2e-10*n);
    for k=[-3 0 1 3 n]
        verifyLessThan(testCase,abs(S.finite.moments(k)-trace(U^k)/n),2e-10*n);
    end
end
end

function test_rank_one_errors(testCase)
verifyError(testCase,@()unitary_rank_one_update(4,.3,ones(3,1)), ...
    'unitary_rank_one_update:Alpha');
verifyError(testCase,@()unitary_rank_one_update(4,.3,zeros(3,1),.4), ...
    'unitary_rank_one_update:Beta');
verifyError(testCase,@()unitary_rank_one_update(4,.3,zeros(3,1),1,zeros(4,1)), ...
    'unitary_rank_one_update:Vector');
F=unitary_rank_one_update(4,[],[],[],[],'operator');
verifyError(testCase,@()F(ones(3,1),'notransp'),'um_rank_one_apply:Dimension');
verifyError(testCase,@()F(ones(4,1),'bad'),'um_rank_one_apply:Flag');
end

function test_etna_all_types_are_unitary_hessenberg(testCase)
for kind={'I','II','III'}
    name=kind{1};
    for n=[8 12 16]
        U=etna_hessenberg(n,name,11);
        verifyEqual(testCase,size(U),[n n]);
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),3e-11*n);
        verifyLessThan(testCase,norm(tril(U,-2),'fro'),3e-11*n);
        V=etna_hessenberg(n,name,11);
        verifyLessThan(testCase,norm(U-V,'fro'),1e-14);
        W=etna_hessenberg(n,name,12);
        verifyGreaterThan(testCase,norm(U-W,'fro'),1e-6);
        if strcmp(name,'II'),verifyTrue(testCase,isreal(U));end
        S=etna_hessenberg(n,name,11,[],'spectral');
        verifyEqual(testCase,S.finite.kind,'unknown');
        F=etna_hessenberg(n,name,11,[],'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(testCase,F([],'size'),[n n]);
        verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),2e-11*n);
        verifyLessThan(testCase,norm(F(X,'transp')-U'*X,'fro'),2e-11*n);
        verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-11*n);
    end
end
end

function test_etna_type_three_coupling_and_block_limit(testCase)
n=24;B=n/4;
small=1e-6;
U=etna_hessenberg(n,'III',7,small);
larger=etna_hessenberg(n,'III',7,.02);
for k=[B 2*B 3*B]
    verifyLessThan(testCase,abs(abs(U(k+1,k))-small),1e-10);
    verifyLessThan(testCase,abs(abs(larger(k+1,k))-.02),1e-11);
end
ref=U(1:B,1:B);
for j=2:4
    idx=(j-1)*B+(1:B);
    verifyLessThan(testCase,norm(U(idx,idx)-ref,'fro'),2e-6);
end
verifyLessThan(testCase,norm(U-kron(eye(4),ref),'fro'),3e-6);
end

function test_etna_rng_preservation_and_errors(testCase)
old=rng;
etna_hessenberg(16,'I',21);
verifyEqual(testCase,rng,old);
etna_hessenberg(16,'II',3,[],'spectral');
verifyEqual(testCase,rng,old);
verifyError(testCase,@()etna_hessenberg(10,'III',1),'etna_hessenberg:Order');
verifyError(testCase,@()etna_hessenberg(8,'IV',1),'etna_hessenberg:Kind');
verifyError(testCase,@()etna_hessenberg(8,'III',1,0),'etna_hessenberg:Coupling');
verifyError(testCase,@()etna_hessenberg(8,'III',1,1),'etna_hessenberg:Coupling');
verifyError(testCase,@()etna_hessenberg(8,'II',1,[],'bad'),'etna_hessenberg:Mode');
end
