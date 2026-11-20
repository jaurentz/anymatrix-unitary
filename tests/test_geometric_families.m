function tests=test_geometric_families
%TEST_GEOMETRIC_FAMILIES Independent geometric and spectral regression tests.
tests=functiontests(localfunctions);
end

function test_sphere_parallel_transport_and_actions(testCase)
for n=[2 3 9]
    theta=.69;p=zeros(n,1);p(1)=1;
    q=zeros(n,1);q(1)=cos(theta);q(2)=sin(theta);
    R=sphere_geodesic(n,p,q);
    verifyLessThan(testCase,norm(R'*R-eye(n),'fro'),3e-12*n);
    verifyLessThan(testCase,norm(R*p-q),3e-12*n);
    w=zeros(n,1);w(2)=1;
    transported=w-((w'*q)/(1+p'*q))*(p+q);
    verifyLessThan(testCase,norm(R*w-transported),3e-12*n);
    S=sphere_geodesic(n,p,q,'spectral');
    verifyEqual(testCase,S.finite.kind,'exact');
    assert_spectrum(testCase,R,S.finite.eigenvalues());
    F=sphere_geodesic(n,p,q,'operator');
    verify_operator(testCase,R,F);
end
end

function test_grassmann_geodesic_and_complex_frames(testCase)
n=11;angles=[.2;.72;1.4];B=eye(n,6);
G=grassmann_geodesic(n,angles,B);
verifyLessThan(testCase,norm(G'*G-eye(n),'fro'),3e-12*n);
expected=B(:,1:3)*diag(cos(angles))+B(:,4:6)*diag(sin(angles));
verifyLessThan(testCase,norm(G*B(:,1:3)-expected,'fro'),3e-12*n);
S=grassmann_geodesic(n,angles,B,'spectral');
assert_spectrum(testCase,G,S.finite.eigenvalues());
verify_operator(testCase,G,grassmann_geodesic(n,angles,B,'operator'));
[Q,~]=qr(randn(n,6)+1i*randn(n,6),0);
Gc=grassmann_geodesic(n,angles,Q);
verifyLessThan(testCase,norm(Gc'*Gc-eye(n),'fro'),3e-12*n);
verify_operator(testCase,Gc,grassmann_geodesic(n,angles,Q,'operator'));
end

function test_geodesic_triangle_holonomy(testCase)
for n=[3 8]
    P=eye(n,3);
    H=geodesic_holonomy(n,P);
    verifyLessThan(testCase,norm(H'*H-eye(n),'fro'),2e-12*n);
    verifyLessThan(testCase,norm(H*P(:,1)-P(:,1)),2e-12*n);
    S=geodesic_holonomy(n,P,'spectral');
    verifyEqual(testCase,S.finite.kind,'exact');
    assert_spectrum(testCase,H,S.finite.eigenvalues());
    verifyLessThan(testCase,abs(S.finite.moments(1)-(n-2)/n),2e-12);
    verify_operator(testCase,H,geodesic_holonomy(n,P,'operator'));
end
P=eye(5,4);H=geodesic_holonomy(5,P);
S=geodesic_holonomy(5,P,'spectral');
verifyEqual(testCase,S.finite.kind,'unknown');
verifyLessThan(testCase,norm(H'*H-eye(5),'fro'),2e-11);
verify_operator(testCase,H,geodesic_holonomy(5,P,'operator'));
end

function test_connection_cycle_spectrum_sparse_operator(testCase)
for N=[1 2 4 13]
    angles=.13*(1:N)';
    U=connection_cycle(N,angles);
    verifyTrue(testCase,issparse(U));
    verifyLessThan(testCase,norm(U'*U-speye(2*N),'fro'),2e-12*N);
    S=connection_cycle(N,angles,'spectral');
    verifyEqual(testCase,S.dimension,2*N);
    assert_spectrum(testCase,full(U),S.finite.eigenvalues());
    verify_operator(testCase,U,connection_cycle(N,angles,'operator'));
end
N=4096;
F=connection_cycle(N,[],'operator');
X=randn(2*N,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),1e-9);
end

function test_lie_connection_exact_spectrum(testCase)
for n=[1 2 4 7]
    kappa=.17*(0:n-1)'.^2;t=1.4;
    Q=lie_geodesic_transport(n,kappa,t);
    verifyTrue(testCase,issparse(Q));
    verifyLessThan(testCase,norm(Q'*Q-speye(n*n),'fro'),2e-11*n);
    S=lie_geodesic_transport(n,kappa,t,'spectral');
    verifyEqual(testCase,S.dimension,n*n);
    assert_spectrum(testCase,full(Q),S.finite.eigenvalues());
    for k=[0 1 2 3]
        verifyLessThan(testCase,abs(S.finite.moments(k)-trace(Q^k)/(n*n)),2e-11*n);
    end
    verify_operator(testCase,Q,lie_geodesic_transport(n,kappa,t,'operator'));
end
end

function test_geometric_validation(testCase)
p=[1;0;0];
verifyError(testCase,@()sphere_geodesic(3,p,-p),'sphere_geodesic:Antipodal');
verifyError(testCase,@()grassmann_geodesic(5,[.1;.2],ones(5,4)), ...
    'grassmann_geodesic:Frames');
verifyError(testCase,@()connection_cycle(4,[1;2;3;1i]), ...
    'connection_cycle:Angles');
verifyError(testCase,@()geodesic_holonomy(3,[p -p [0;1;0]]), ...
    'geodesic_holonomy:Antipodal');
verifyError(testCase,@()lie_geodesic_transport(3,[1 2]), ...
    'lie_geodesic_transport:Kappa');
end

function verify_operator(testCase,A,F)
n=size(A,1);X=randn(n,3)+1i*randn(n,3);
verifyEqual(testCase,F([],'size'),[n n]);
verifyLessThan(testCase,norm(F(X,'notransp')-A*X,'fro'),3e-11*n);
verifyLessThan(testCase,norm(F(X,'transp')-A'*X,'fro'),3e-11*n);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),3e-11*n);
end

function assert_spectrum(testCase,A,z)
x=eig(full(A));verifyEqual(testCase,numel(x),numel(z));
verifyLessThan(testCase,max(arrayfun(@(v)min(abs(v-z)),x)),1e-9*numel(z));
verifyLessThan(testCase,max(arrayfun(@(v)min(abs(v-x)),z)),1e-9*numel(z));
end
