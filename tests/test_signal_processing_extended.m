function tests=test_signal_processing_extended
%TEST_SIGNAL_PROCESSING_EXTENDED New FFT and eigensystem transforms.
tests=functiontests(localfunctions);
end

function test_discrete_fourier_exact(testCase)
for n=[1 2 3 4 7 12]
    F=discrete_fourier(n);
    check_unitary(testCase,F,2e-11*n);
    verifyLessThan(testCase,norm(F-F.','fro'),2e-11*n);
    j=0:n-1;
    R=full(sparse(mod(-j,n)+1,j+1,1,n,n));
    verifyLessThan(testCase,norm(F*F-R,'fro'),2e-11*n);
    verifyLessThan(testCase,norm(F^4-eye(n),'fro'),5e-11*n);
    check_operator(testCase,F,discrete_fourier(n,'operator'));
end
end

function test_four_discrete_sine_types(testCase)
for n=[1 2 3 5 8 13]
    S=cell(1,4);
    for kind=1:4
        S{kind}=discrete_sine(n,kind);
        check_unitary(testCase,S{kind},4e-11*n);
        check_operator(testCase,S{kind},discrete_sine(n,kind,'operator'));
        if kind==1 || kind==4
            verifyLessThan(testCase,norm(S{kind}-S{kind}','fro'),3e-11*n);
        end
    end
    verifyLessThan(testCase,norm(S{2}'-S{3},'fro'),4e-11*n);
end
check_operator(testCase,discrete_sine(8),discrete_sine(8,'operator'));
verifyError(testCase,@()discrete_sine(8,5),'MATLAB:expectedLessEqual');
end

function test_walsh_natural_and_sequency(testCase)
for n=[1 2 4 8 16]
    H=walsh_hadamard(n,'natural');
    Q=walsh_hadamard(n,'sequency');
    check_unitary(testCase,H,2e-12*n);
    check_unitary(testCase,Q,2e-12*n);
    verifyLessThan(testCase,norm(H-H','fro'),2e-12*n);
    verifyLessThan(testCase,norm(H^2-eye(n),'fro'),2e-12*n);
    verifyEqual(testCase,unique(abs(H(:))),1/sqrt(n),'AbsTol',2e-12);
    flips=sum(diff(Q,1,2)~=0,2);
    verifyEqual(testCase,flips,(0:n-1)');
    check_operator(testCase,H,walsh_hadamard(n,'natural','operator'));
    check_operator(testCase,Q,walsh_hadamard(n,'operator'));
end
verifyError(testCase,@()walsh_hadamard(5),'walsh_hadamard:Order');
verifyError(testCase,@()walsh_hadamard(8,'bogus'),'walsh_hadamard:Ordering');
end

function test_slepian_concentration_vs_commuting_tridiagonal(testCase)
for n=[1 2 4 8 12 20]
    for W=[.09 .21 .38]
        U=dpss_transform(n,W);
        check_unitary(testCase,U,5e-11*n);
        j=(0:n-1)';
        delta=j-j';
        B=2*W*eye(n);
        I=delta~=0;
        B(I)=sin(2*pi*W*delta(I))./(pi*delta(I));
        C=U*B*U';
        verifyLessThan(testCase,norm(C-diag(diag(C)),'fro'),8e-10*n);
        concentration=diag(C);
        verifyTrue(testCase,all(concentration>=-2e-11));
        verifyTrue(testCase,all(concentration<=1+2e-11));
        verifyTrue(testCase,all(diff(concentration)<=1e-10));
        if n>1
            d=((n-1)/2-j).^2*cos(2*pi*W);
            e=(1:n-1)';e=e.*(n-e)/2;
            T=diag(d)+diag(e,1)+diag(e,-1);
            verifyLessThan(testCase,norm(T*B-B*T,'fro'),2e-9*n);
        end
    end
end
verifyEqual(testCase,dpss_transform(7),dpss_transform(7,.25),'AbsTol',1e-12);
end

function test_graph_laplacian_basis_and_degeneracies(testCase)
for n=[1 2 3 7 12]
    U=graph_fourier(n);
    check_unitary(testCase,U,2e-10*n);
    A=zeros(n);for j=1:n-1,A(j,j+1)=1;A(j+1,j)=1;end
    L=diag(sum(A,2))-A;
    d=sort(eig(L),'ascend');
    verifyLessThan(testCase,norm(U*L*U'-diag(d),'fro'),4e-10*n);
end
for n=[4 6 9]
    A=zeros(n);A(1,2:end)=1;A(2:end,1)=1;
    B=graph_fourier(n,A);
    check_unitary(testCase,B,2e-10*n);
    verifyLessThan(testCase,norm(B-graph_fourier(n,A),'fro'),8e-10*n);
    L=diag(sum(A,2))-A;
    d=sort(eig(L),'ascend');
    verifyLessThan(testCase,norm(B*L*B'-diag(d),'fro'),5e-9*n);
end
A=[0 1;0 0];
verifyError(testCase,@()graph_fourier(2,A),'graph_fourier:Adjacency');
end

function test_implicit_large_transform_actions(testCase)
n=2048;
X=randn(n,2)+1i*randn(n,2);
ops={discrete_fourier(n,'operator'),discrete_sine(n,1,'operator'), ...
    discrete_sine(n,2,'operator'),discrete_sine(n,3,'operator'), ...
    discrete_sine(n,4,'operator'),walsh_hadamard(n,'natural','operator'), ...
    walsh_hadamard(n,'sequency','operator')};
for k=1:numel(ops)
    F=ops{k};
    verifyEqual(testCase,F([],'size'),[n n]);
    verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),3e-9);
end
% DPSS and general graph transforms intentionally remain explicit matrices.
end

function test_operator_bad_shapes(testCase)
F=walsh_hadamard(8,'operator');
verifyError(testCase,@()F(ones(7,2),'notransp'),'um_walsh_apply:Dimension');
verifyError(testCase,@()F(ones(8,2),'bad'),'um_walsh_apply:Flag');
end

function check_unitary(testCase,U,tol)
n=size(U,1);
verifyEqual(testCase,size(U),[n n]);
verifyLessThan(testCase,norm(full(U'*U-eye(n)),'fro'),tol);
end

function check_operator(testCase,A,F)
n=size(A,1);tol=8e-11*max(1,n);
verifyTrue(testCase,isa(F,'function_handle'));
verifyEqual(testCase,F([],'size'),[n n]);
X=randn(n,3)+1i*randn(n,3);
Y=randn(n,3)+1i*randn(n,3);
verifyLessThan(testCase,norm(F(X,'notransp')-A*X,'fro'),tol);
verifyLessThan(testCase,norm(F(Y,'transp')-A'*Y,'fro'),tol);
verifyLessThan(testCase,norm(F(Y,'ctransp')-A'*Y,'fro'),tol);
verifyLessThan(testCase,norm(F(eye(n),'notransp')-A,'fro'),tol);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),tol);
end
