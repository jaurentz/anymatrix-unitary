function tests=test_lohe_matrix_flow
%TEST_LOHE_MATRIX_FLOW Nonlinear unitary flow, invariances and order.
tests=functiontests(localfunctions);
end

function test_inventory_and_default_rng_reproducibility(testCase)
[empty,props]=lohe_matrix_flow();
verifyEmpty(testCase,empty);
verifyTrue(testCase,any(strcmp(props,'unitary')));
state=rng;
[A,~,All]=lohe_matrix_flow(3);
verifyEqual(testCase,rng,state);
[B,~,AllB]=lohe_matrix_flow(3);
verifyEqual(testCase,A,B,'AbsTol',0);
verifyEqual(testCase,All,AllB,'AbsTol',0);
verifyEqual(testCase,size(All),[3 3 3]);
verifyEqual(testCase,A,All(:,:,1));
for j=1:3
    verifyLessThan(testCase,norm(All(:,:,j)'*All(:,:,j)-eye(3),'fro'),4e-13);
end
end

function test_exact_decoupled_and_single_oscillator(testCase)
d=3;n=4;
H=complex(zeros(d,d,n));
U0=complex(zeros(d,d,n));
for j=1:n
    T=randn(d)+1i*randn(d);
    H(:,:,j)=(T+T')/2;
    U0(:,:,j)=haar_unitary(d,j+25);
end
time=.7;
for steps=[1 2 11]
    [A,~,ensemble]=lohe_matrix_flow(d,n,0,time,H,U0,2,steps);
    verifyEqual(testCase,A,ensemble(:,:,2),'AbsTol',5e-14);
    for j=1:n
        expected=expm(-1i*time*H(:,:,j))*U0(:,:,j);
        verifyLessThan(testCase,norm(ensemble(:,:,j)-expected,'fro'), ...
            2e-13);
    end
end
V=lohe_matrix_flow(d,1,4,time,H(:,:,1),U0(:,:,1),1,5);
verifyLessThan(testCase,norm(V-expm(-1i*time*H(:,:,1))*U0(:,:,1),'fro'), ...
    2e-13);
end

function test_synchronized_exact_propagation(testCase)
d=3;n=5;
H1=[1 .2+.1i 0;.2-.1i 0 -.4i;0 .4i -1];
V=haar_unitary(d,57);
H=repmat(H1,1,1,n);
U0=repmat(V,1,1,n);
[U,~,All]=lohe_matrix_flow(d,n,3,0.8,H,U0,3,8);
exact=expm(-.8i*H1)*V;
verifyLessThan(testCase,norm(U-exact,'fro'),5e-13);
for j=1:n
    verifyLessThan(testCase,norm(All(:,:,j)-exact,'fro'),5e-13);
end
end

function test_right_unitary_covariance(testCase)
d=3;n=3;
U0=complex(zeros(d,d,n));
H=complex(zeros(d,d,n));
for j=1:n
    U0(:,:,j)=haar_unitary(d,50+j);
    T=randn(d)+1i*randn(d);
    H(:,:,j)=(T+T')/2;
end
R=haar_unitary(d,808);
Urot=complex(zeros(d,d,n));
for j=1:n,Urot(:,:,j)=U0(:,:,j)*R;end
A=lohe_matrix_flow(d,n,.65,.7,H,U0,2,25,'ensemble');
B=lohe_matrix_flow(d,n,.65,.7,H,Urot,2,25,'ensemble');
for j=1:n
    verifyLessThan(testCase,norm(B(:,:,j)-A(:,:,j)*R,'fro'),2e-12);
end
end

function test_scalar_Kuramoto_exact_two_oscillator(testCase)
% With H=0, theta_1=-a, theta_2=a:
% theta_2(t)=atan(exp(-kappa*t)*tan(a)), theta_1=-theta_2.
a=.6; time=1.1;kappa=1.2;
U0=reshape(exp(1i*[-a a]),[1,1,2]);
H=zeros(1,1,2);
beta=atan(exp(-kappa*time)*tan(a));
exact=exp(1i*[-beta beta]);
errors=zeros(1,3);
steps=[12 24 48];
for j=1:3
    U=lohe_matrix_flow(1,2,kappa,time,H,U0,1,steps(j),'ensemble');
    errors(j)=norm(reshape(U,1,2)-exact);
    verifyLessThan(testCase,max(abs(abs(U(:))-1)),1e-13);
end
verifyLessThan(testCase,errors(3),3e-4);
verifyLessThan(testCase,errors(2),errors(1)/3);
verifyLessThan(testCase,errors(3),errors(2)/3);
% SU(1) is not preserved when coupled: U_1(0)=1 evolves away from 1.
U0=reshape([1 exp(1i*.8)],[1,1,2]);
U=lohe_matrix_flow(1,2,1,.4,H,U0,1,30);
verifyGreaterThan(testCase,abs(U-1),.02);
end

function test_real_orthogonal_subfamily(testCase)
d=2;n=3;
theta=[-.3 .5 1.2];
U0=zeros(d,d,n);
for j=1:n
    c=cos(theta(j));s=sin(theta(j));
    U0(:,:,j)=[c -s;s c];
end
H=zeros(d,d,n);
U=lohe_matrix_flow(d,n,2,1,H,U0,1,24,'ensemble');
for j=1:n
    X=U(:,:,j);
    verifyLessThan(testCase,norm(imag(X),'fro'),1e-13);
    verifyLessThan(testCase,norm(X.'*X-eye(d),'fro'),2e-13);
    verifyLessThan(testCase,abs(det(X)-1),2e-13);
end
end

function test_initial_state_and_unknown_spectrum(testCase)
d=3;
S=lohe_matrix_flow(d,[],[],[],[],[],[],[],'spectral');
verifyEqual(testCase,S.dimension,d);
verifyEqual(testCase,S.finite.kind,'unknown');
verifyEqual(testCase,S.ensemble.kind,'unknown');
verifyEqual(testCase,S.family,'lohe_matrix_flow');
A=lohe_matrix_flow(d,3,5,0,[],[],2,8);
B=lohe_matrix_flow(d,3,0,0,[],[],2,8);
verifyEqual(testCase,A,B,'AbsTol',0);
end

function test_input_validation(testCase)
verifyError(testCase,@()lohe_matrix_flow(2,3,-1), ...
    'MATLAB:expectedNonnegative');
verifyError(testCase,@()lohe_matrix_flow(2,2,1,.5,cat(3,[1 2;0 1],eye(2))), ...
    'lohe_matrix_flow:Hermitian');
bad=cat(3,eye(2),ones(2));
verifyError(testCase,@()lohe_matrix_flow(2,2,1,.5,[],bad), ...
    'lohe_matrix_flow:Unitary');
verifyError(testCase,@()lohe_matrix_flow(2,2,1,.5,[],[],3), ...
    'MATLAB:notLessEqual');
verifyError(testCase,@()lohe_matrix_flow(2,2,1,.5,[],[],1,0), ...
    'MATLAB:expectedPositive');
verifyError(testCase,@()lohe_matrix_flow(2,2,1,.5,[],[],1,8,'operator'), ...
    'lohe_matrix_flow:Mode');
end
