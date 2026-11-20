function tests=test_geodesic_extensions
%TEST_GEODESIC_EXTENSIONS Stiefel, Grassmann holonomy and Wigner tests.
tests=functiontests(localfunctions);
end

function test_stiefel_geodesic_tangent_and_exact_reduction(t)
for n=[6 12]
    p=2;A=[0 -.4;.4 0];R=diag([.27 .65]);B=eye(n,2*p);
    G=stiefel_geodesic(n,p,A,R,B,1);
    verifyLessThan(t,norm(G'*G-eye(n),'fro'),2e-11*n);
    tangent=B(:,1:p)*A+B(:,p+1:2*p)*R;
    h=1e-5;
    Gp=stiefel_geodesic(n,p,A,R,B,h);
    Gm=stiefel_geodesic(n,p,A,R,B,-h);
    verifyLessThan(t,norm((Gp*B(:,1:p)-Gm*B(:,1:p))/(2*h)-tangent,'fro'),2e-9);
    S=stiefel_geodesic(n,p,A,R,B,1,'spectral');
    verifyEqual(t,S.finite.kind,'exact');
    assert_spec(t,G,S.finite.eigenvalues());
    F=stiefel_geodesic(n,p,A,R,B,1,'operator');
    verify_op(t,G,F);
end
n=9;p=2;[B,~]=qr(randn(n,2*p)+1i*randn(n,2*p),0);
A=1i*diag([.3 -.5]);R=[.2 .4i;.7 .1];
G=stiefel_geodesic(n,p,A,R,B,.7);
verifyLessThan(t,norm(G'*G-eye(n),'fro'),3e-11*n);
assert_spec(t,G,stiefel_geodesic(n,p,A,R,B,.7,'spectral').finite.eigenvalues());
verify_op(t,G,stiefel_geodesic(n,p,A,R,B,.7,'operator'));
end

function test_grassmann_nonabelian_gauge_covariance(t)
n=5;r=2;a=.51;b=.63;
P=zeros(n,r,3);
P(:,:,1)=eye(n,r);
P(:,:,2)=P(:,:,1);P(1,1,2)=cos(a);P(3,1,2)=sin(a);
P(:,:,3)=P(:,:,1);P(2,2,3)=cos(b);P(3,2,3)=exp(.7i)*sin(b);
H=grassmann_holonomy(n,r,P);
verifyLessThan(t,norm(H'*H-eye(r),'fro'),2e-12);
verifyGreaterThan(t,norm(H-eye(r),'fro'),1e-4);
S=grassmann_holonomy(n,r,P,'spectral');
verifyEqual(t,S.finite.kind,'unknown');
Q=zeros(r,r,3);P2=P;
for k=1:3
    [Q(:,:,k),~]=qr(randn(r)+1i*randn(r));
    P2(:,:,k)=P(:,:,k)*Q(:,:,k);
end
H2=grassmann_holonomy(n,r,P2);
verifyLessThan(t,norm(H2-Q(:,:,1)'*H*Q(:,:,1),'fro'),2e-11);
assert_spec(t,H,H2);
end

function test_grassmann_rank_one_bargmann(t)
p=[1;0;0];q=[cos(.4);sin(.4);0];
v=[cos(.5);0;exp(.8i)*sin(.5)];
P=cat(3,p,q,v);
H=grassmann_holonomy(3,1,P);
phase=(q'*p)*(v'*q)*(p'*v);
phase=phase/abs(phase);
verifyLessThan(t,abs(H-phase),2e-12);
S=grassmann_holonomy(3,1,P,'spectral');
verifyEqual(t,S.finite.kind,'exact');
verifyLessThan(t,abs(S.finite.eigenvalues()-H),2e-12);
verifyLessThan(t,abs(S.finite.moments(5)-H^5),2e-12);
end

function test_wigner_representation_and_spectrum(t)
for n=[1 2 3 5 8]
    a=.37;b=.91;c=-1.2;
    U=wigner_geodesic(n,a,b,c);
    verifyLessThan(t,norm(U'*U-eye(n),'fro'),4e-11*n);
    S=wigner_geodesic(n,a,b,c,'spectral');
    verifyEqual(t,S.finite.kind,'exact');
    assert_spec(t,U,S.finite.eigenvalues());
    verifyLessThan(t,abs(det(U)-1),3e-10*n);
    for k=[0 1 2 -2]
        verifyLessThan(t,abs(S.finite.moments(k)-trace(U^k)/n),2e-10*n);
    end
    D=wigner_geodesic(n,a,0,c);
    m=(-(n-1)/2:(n-1)/2)';
    verifyLessThan(t,norm(D-diag(exp(-1i*(a+c)*m)),'fro'),2e-12*n);
end
% Fundamental SU(2) representation must retain a nontrivial 2pi lift.
U=wigner_geodesic(2,2*pi,0,0);
verifyLessThan(t,norm(U+eye(2),'fro'),2e-12);
end

function test_rejections(t)
verifyError(t,@()stiefel_geodesic(4,3),'stiefel_geodesic:Order');
verifyError(t,@()stiefel_geodesic(6,2,eye(2)),'stiefel_geodesic:Skew');
verifyError(t,@()grassmann_holonomy(3,3),'grassmann_holonomy:Order');
P=cat(3,[1;0;0],[0;1;0],[0;0;1]);
verifyError(t,@()grassmann_holonomy(3,1,P),'grassmann_holonomy:CutLocus');
verifyError(t,@()wigner_geodesic(4,.1,.2,.3,'operator'),'wigner_geodesic:Mode');
end

function assert_spec(t,A,z)
v=eig(full(A));z=z(:);
verifyEqual(t,numel(v),numel(z));
verifyLessThan(t,max(arrayfun(@(x)min(abs(x-z)),v)),2e-9*numel(z));
verifyLessThan(t,max(arrayfun(@(x)min(abs(x-v)),z)),2e-9*numel(z));
end

function verify_op(t,A,F)
n=size(A,1);X=randn(n,3)+1i*randn(n,3);
verifyEqual(t,F([],'size'),[n n]);
verifyLessThan(t,norm(F(X,'notransp')-A*X,'fro'),3e-11*n);
verifyLessThan(t,norm(F(X,'transp')-A'*X,'fro'),3e-11*n);
verifyLessThan(t,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),3e-11*n);
end
