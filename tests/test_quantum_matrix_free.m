function tests=test_quantum_matrix_free
%TEST_QUANTUM_MATRIX_FREE New bosonic, cavity-QED and quantum-network actions.
%   Also regression-tests the shared tridiagonal exponential on kicked top.
tests=functiontests(localfunctions);
end

function test_tridiagonal_complex_exponential(testCase)
s=rng;cleanup=onCleanup(@()rng(s)); %#ok<NASGU>
rng(11,'twister');
for n=[1 2 3 7 17 70]
    d=.35*sin((1:n)'*.31);
    u=.2*exp(1i*(1:n-1)'*.41);
    H=diag(d)+diag(u,1)+diag(conj(u),-1);
    X=randn(n,3)+1i*randn(n,3);
    for time=[0 .13 -.86 2.5]
        result=um_tridiag_exp_apply(X,d,u,time);
        reference=expm(-1i*time*H)*X;
        verifyLessThan(testCase,norm(result-reference,'fro'), ...
            3e-11*max(1,norm(X,'fro')));
    end
end
end

function test_bose_matrix_free_vs_explicit(testCase)
for n=[0 1 2 3 8 16]
    for parameters={[.7 .4 .81 1.27 .3],[1.1 .05 -.6 .5 -.2], ...
                   [0 .2 1 .7 .9],[.6 0 0 0 0]}
        p=parameters{1};
        A=kicked_bose_hubbard(n,p(1),p(2),p(3),p(4),p(5));
        O=kicked_bose_hubbard(n,p(1),p(2),p(3),p(4),p(5),'operator');
        check_operator(testCase,A,O);
        if n<=2
            P=kicked_bose_hubbard(n,p(1),p(2),p(3),p(4),p(5),'implicit');
            check_operator(testCase,A,P);
        end
    end
end
end

function test_kicked_top_shared_spin_action(testCase)
for n=[2 3 5 11 31]
    for twist=[0 .35 2.7]
        for p=[0 -.8 .76]
            A=kicked_top(n,twist,p);
            O=kicked_top(n,twist,p,'operator');
            check_operator(testCase,A,O);
        end
    end
end
end

function test_jaynes_cummings_blocks_and_edges(testCase)
for n=[0 1 2 7 31]
    for g=[0 .3 1.2]
        for delta=[-.4 0 .7]
            A=jaynes_cummings(n,g,delta,1.11,1.7);
            verifyTrue(testCase,issparse(A));
            O=jaynes_cummings(n,g,delta,1.11,1.7,'operator');
            check_operator(testCase,A,O);
        end
    end
end
end

function test_chalker_coddington_network_vs_sparse(testCase)
for dims=[1 1;1 3;3 1;2 2;3 4]'
    nx=dims(1);ny=dims(2);
    n=4*nx*ny;
    phi=.17*(1:n)'+.35*sin((1:n)');
    for angles={0,pi/4,[.2 .87],pi/2}
        A=chalker_coddington(nx,ny,angles{1},phi);
        O=chalker_coddington(nx,ny,angles{1},phi,'operator');
        check_operator(testCase,A,O);
    end
end
end

function test_disjoint_cc_link_ports(testCase)
for nx=[1 2 3]
    for ny=[1 2 4]
        [a,b,gates]=um_cc_geometry(nx,ny,[.3;.71]);
        n=4*nx*ny;
        verifyEqual(testCase,size(a),[2 n/2]);
        verifyEqual(testCase,size(b),[2 n/2]);
        verifyEqual(testCase,sort(a(:)),(1:n)');
        verifyEqual(testCase,sort(b(:)),(1:n)');
        verifyEqual(testCase,size(gates),[4 n/2]);
        A=chalker_coddington(nx,ny,[.3 .71]);
        check_operator(testCase,A,chalker_coddington(nx,ny,[.3 .71],0,'operator'));
    end
end
end

function test_large_implicit_actions(testCase)
s=rng;cleanup=onCleanup(@()rng(s)); %#ok<NASGU>
rng(19,'twister');
% A dense matrix is not constructed for these cases.
B=kicked_bose_hubbard(256,.7,.2,.51,1,.16,'operator');
check_adjoint(testCase,B,257,2e-10);
T=kicked_top(257,2,.48,'operator');
check_adjoint(testCase,T,257,2e-10);
J=jaynes_cummings(4096,.8,.31,.37,1.5,'operator');
check_adjoint(testCase,J,8194,3e-11);
C=chalker_coddington(50,30,[.3 .7],.27,'operator');
check_adjoint(testCase,C,6000,3e-11);
end

function test_action_flags_and_dimensions(testCase)
B=kicked_bose_hubbard(3,[],[],[],[],[],'operator');
verifyError(testCase,@()B(ones(5,1),'notransp'), ...
    'um_kick_tridiag_apply:Dimension');
verifyError(testCase,@()B(eye(4),'invalid'), ...
    'um_kick_tridiag_apply:Flag');
T=kicked_top(4,[],[],'operator');
verifyError(testCase,@()T(ones(3,1),'notransp'), ...
    'um_kick_tridiag_apply:Dimension');
J=jaynes_cummings(3,[],[],[],[],'operator');
verifyError(testCase,@()J(ones(9,1),'notransp'), ...
    'um_jaynes_apply:Dimension');
verifyError(testCase,@()J(eye(8),'invalid'), ...
    'um_jaynes_apply:Flag');
C=chalker_coddington(2,3,[],[],'operator');
verifyError(testCase,@()C(ones(25,1),'notransp'), ...
    'um_scattering_apply:Dimension');
verifyError(testCase,@()C(eye(24),'invalid'), ...
    'um_scattering_apply:Flag');
verifyError(testCase,@()jaynes_cummings(3,[],[],[],[],'bogus'), ...
    'um_is_operator:Mode');
end

function check_operator(testCase,A,O)
n=size(A,1);
verifyTrue(testCase,isa(O,'function_handle'));
verifyEqual(testCase,O([],'size'),[n n]);
X=reshape(sin((1:n*3)*.71),n,3)+1i*reshape(cos((1:n*3)*.39),n,3);
Z=reshape(cos((1:n*3)*.43),n,3)+1i*reshape(sin((1:n*3)*.29),n,3);
tol=5e-11*max(1,n);
verifyLessThan(testCase,norm(O(X,'notransp')-A*X,'fro'),tol);
verifyLessThan(testCase,norm(O(Z,'transp')-A'*Z,'fro'),tol);
verifyLessThan(testCase,norm(O(Z,'ctransp')-A'*Z,'fro'),tol);
verifyLessThan(testCase,norm(O(eye(n),'notransp')-full(A),'fro'),tol);
verifyLessThan(testCase,norm(O(O(X,'notransp'),'transp')-X,'fro'),tol);
end

function check_adjoint(testCase,O,n,tol)
verifyEqual(testCase,O([],'size'),[n n]);
X=randn(n,2)+1i*randn(n,2);
Y=randn(n,2)+1i*randn(n,2);
Z=O(X,'notransp');
W=O(Y,'transp');
verifyLessThan(testCase,norm(O(Z,'transp')-X,'fro'),tol*norm(X,'fro'));
verifyLessThan(testCase,abs(sum(conj(Z(:)).*Y(:))-sum(conj(X(:)).*W(:))), ...
    tol*norm(X,'fro')*norm(Y,'fro'));
end
