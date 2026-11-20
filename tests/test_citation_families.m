function tests=test_citation_families
%TEST_CITATION_FAMILIES Citation-led unitary construction regressions.
tests=functiontests(localfunctions);
end

function test_rational_opuc_is_unitary(testCase)
for n=[1 2 5 12]
    beta=exp(1i*.37);
    schur=(.08+.12i)*ones(n-1,1);
    poles=(1.7+.2i)*ones(n-1,1);
    U=rational_opuc(n,poles,schur,beta);
    verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),3e-11*n);
    verifyLessThan(testCase,max(abs(abs(eig(U))-1)),3e-11);
    H=opuc_matrix(n,schur,beta,'hessenberg');
    verifyLessThan(testCase,norm(rational_opuc(n,inf(n-1,1),schur,beta)-H,'fro'),2e-11*n);
end
end

function test_rational_opuc_noncommuting_poles(testCase)
n=9;
poles=[inf;2+1i;-1.2+.4i;2i;1.7;inf;-2i;3+.1i];
schur=.25*exp(1i*(1:n-1)'*.31);
U=rational_opuc(n,poles,schur,exp(1i*.47));
verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),3e-11*n);
S=rational_opuc(n,poles,schur,exp(1i*.47),'spectral');
verifyEqual(testCase,S.dimension,n);
verifyEqual(testCase,S.finite.kind,'unknown');
verifyFalse(testCase,isempty(S.reference));
end

function test_rational_opuc_input_validation(testCase)
verifyError(testCase,@()rational_opuc(4,[2;0;inf]),'rational_opuc:Poles');
verifyError(testCase,@()rational_opuc(4,[2;1;inf]),'rational_opuc:Poles');
verifyError(testCase,@()rational_opuc(4,[2;3]),'rational_opuc:Poles');
verifyError(testCase,@()rational_opuc(3,[2;inf],[0;1]),'rational_opuc:Schur');
verifyError(testCase,@()rational_opuc(3,[2;inf],[0;0],.9),'rational_opuc:Boundary');
verifyError(testCase,@()rational_opuc(3,[],[],[],'operator'),'rational_opuc:Mode');
end

function test_blaschke_unitarity_and_determinant(testCase)
for n=[1 2 5 8]
    poles=[.15+.1i;-.35+.2i;.5-.2i];
    V=complex(zeros(n,numel(poles)));
    for k=1:numel(poles)
        V(mod(k-1,n)+1,k)=1;
        V(mod(k,n)+1,k)=V(mod(k,n)+1,k)+.7i;
    end
    % Use a genuine N-by-N unitary factor as the initial matrix.
    [Q,~]=qr(randn(n)+1i*randn(n));
    for omega=[-.8 0 .9 pi]
        z=exp(1i*omega);
        b=(z-poles)./(1-conj(poles)*z);
        U=blaschke_potapov(n,poles,V,omega,Q);
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),5e-11*n);
        verifyLessThan(testCase,abs(det(U)-det(Q)*prod(b)),5e-10*n);
        F=blaschke_potapov(n,poles,V,omega,Q,'operator');
        verifyEqual(testCase,F([],'size'),[n n]);
        X=randn(n,3)+1i*randn(n,3);
        verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),3e-11*n);
        verifyLessThan(testCase,norm(F(X,'transp')-U'*X,'fro'),3e-11*n);
        verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),3e-11*n);
    end
end
end

function test_blaschke_boundaries_and_errors(testCase)
U=blaschke_potapov(6);
verifyLessThan(testCase,norm(U'*U-eye(6),'fro'),5e-11);
S=blaschke_potapov(6,[],[],[],[],'spectral');
verifyEqual(testCase,S.dimension,6);
verifyEqual(testCase,S.finite.kind,'unknown');
verifyError(testCase,@()blaschke_potapov(4,1),'blaschke_potapov:Poles');
verifyError(testCase,@()blaschke_potapov(4,.1,zeros(4,1)), ...
    'blaschke_potapov:Directions');
verifyError(testCase,@()blaschke_potapov(4,.1,ones(3,1)), ...
    'blaschke_potapov:Directions');
verifyError(testCase,@()blaschke_potapov(4,.1,[],0,zeros(4)), ...
    'blaschke_potapov:U0');
F=blaschke_potapov(4,[],[],[],[],'operator');
verifyError(testCase,@()F(ones(3,1),'notransp'),'um_blaschke_apply:Dimension');
verifyError(testCase,@()F(ones(4,1),'invalid'),'um_blaschke_apply:Flag');
end

function test_butterfly_recursive_orthogonality(testCase)
for n=[1 2 4 8 16 32]
    angles=.14*(1:n-1)';
    Q=orthogonal_butterfly(n,angles);
    verifyTrue(testCase,isreal(Q));
    verifyLessThan(testCase,norm(Q'*Q-eye(n),'fro'),5e-12*n);
    verifyLessThan(testCase,abs(det(Q)-1),2e-10*n);
    F=orthogonal_butterfly(n,angles,'operator');
    verifyEqual(testCase,F([],'size'),[n n]);
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(testCase,norm(F(X,'notransp')-Q*X,'fro'),5e-12*n);
    verifyLessThan(testCase,norm(F(X,'transp')-Q'*X,'fro'),5e-12*n);
    verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),5e-12*n);
end
theta=.37;
Q=orthogonal_butterfly(2,theta);
verifyLessThan(testCase,norm(Q-[cos(theta) sin(theta); ...
    -sin(theta) cos(theta)],'fro'),1e-14);
end

function test_butterfly_large_matrix_free(testCase)
n=8192;
angles=.23*ones(n-1,1);
F=orthogonal_butterfly(n,angles,'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),1e-9);
S=orthogonal_butterfly(n,angles,'spectral');
verifyEqual(testCase,S.dimension,n);
verifyEqual(testCase,S.finite.kind,'unknown');
end

function test_butterfly_validation(testCase)
verifyError(testCase,@()orthogonal_butterfly(7),'orthogonal_butterfly:Order');
verifyError(testCase,@()orthogonal_butterfly(4,[.2 .4]), ...
    'orthogonal_butterfly:Angles');
verifyError(testCase,@()orthogonal_butterfly(4,[0 1i 0]), ...
    'orthogonal_butterfly:Angles');
F=orthogonal_butterfly(4,[],'operator');
verifyError(testCase,@()F(ones(3,1),'notransp'),'um_butterfly_apply:Dimension');
verifyError(testCase,@()F(ones(4,1),'invalid'),'um_butterfly_apply:Flag');
end
