function tests=test_geometry_holonomy_wilson
%TEST_GEOMETRY_HOLONOMY_WILSON Closed-loop geometry and gauge invariants.
tests=functiontests(localfunctions);
end

function test_sphere_octant_triangle_holonomy(testCase)
V=eye(3);
H=sphere_geodesic_holonomy(3,V);
expected=[1 0 0;0 0 -1;0 1 0];
verifyLessThan(testCase,norm(H-expected,'fro'),3e-14);
verifyLessThan(testCase,norm(H'*H-eye(3),'fro'),3e-14);
verifyLessThan(testCase,abs(det(H)-1),3e-14);
S=sphere_geodesic_holonomy(3,V,'spectral');
verifyEqual(testCase,S.finite.kind,'exact');
z=S.finite.eigenvalues();
verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),eig(H))),3e-14);
verifyLessThan(testCase,abs(S.finite.moments(1)-1/3),2e-14);
verifyLessThan(testCase,abs(S.finite.moments(2)+1/3),2e-14);
% Reversing the oriented geodesic triangle reverses the SO(2) holonomy.
B=sphere_geodesic_holonomy(3,V(:,[1 3 2]));
verifyLessThan(testCase,norm(B-H','fro'),3e-14);
end

function test_sphere_triangle_holonomy_in_higher_dimensions(testCase)
for n=[3 4 7 11]
    V=zeros(n,3);
    V(1,1)=1;
    V(1,2)=.3;V(2,2)=sqrt(1-.3^2);
    V(1,3)=.2;V(2,3)=.4;V(3,3)=sqrt(1-.2^2-.4^2);
    H=sphere_geodesic_holonomy(n,V);
    verifyLessThan(testCase,norm(H'*H-eye(n),'fro'),5e-12*n);
    verifyLessThan(testCase,norm(H*V(:,1)-V(:,1)),5e-12*n);
    verifyLessThan(testCase,abs(det(H)-1),5e-11*n);
    verifyLessThan(testCase(norm(H-eye(n),'fro')>0),1);
    verifyLessThan(testCase,svd(H-eye(n),'econ')(3),1e-11*n);
    S=sphere_geodesic_holonomy(n,V,'spectral');
    verifyEqual(testCase,S.finite.kind,'exact');
    verifyEqual(testCase,numel(S.finite.eigenvalues()),n);
    z=S.finite.eigenvalues();
    e=eig(H);
    verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),e)),5e-11*n);
    for k=[-5 -2 0 1 2 4]
        verifyLessThan(testCase,abs(S.finite.moments(k)-trace(H^k)/n),5e-11*n);
    end
    verifyLessThan(testCase,abs(sum(S.finite.atoms.weights)-1),3e-15);
end
end

function test_sphere_holonomy_polygon_nonabelian_and_inverse(testCase)
n=5;V=eye(5,4);
H=sphere_geodesic_holonomy(n,V);
verifyLessThan(testCase,norm(H'*H-eye(n),'fro'),2e-11*n);
verifyLessThan(testCase,norm(H*V(:,1)-V(:,1)),2e-11*n);
S=sphere_geodesic_holonomy(n,V,'spectral');
verifyEqual(testCase,S.finite.kind,'unknown');
R=sphere_geodesic_holonomy(n,V(:,[1 4 3 2]));
verifyLessThan(testCase,norm(R-H','fro'),2e-11*n);
F=sphere_geodesic_holonomy(n,V,'operator');
verifyEqual(testCase,F([],'size'),[n n]);
X=randn(n,3)+1i*randn(n,3);
verifyLessThan(testCase,norm(F(X,'notransp')-H*X,'fro'),2e-11*n);
verifyLessThan(testCase,norm(F(X,'ctransp')-H'*X,'fro'),2e-11*n);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-11*n);
end

function test_sphere_holonomy_large_operator_and_validation(testCase)
n=4096;
F=sphere_geodesic_holonomy(n,[],'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-9);
verifyError(testCase,@()sphere_geodesic_holonomy(2),'MATLAB:notGreaterEqual');
verifyError(testCase,@()sphere_geodesic_holonomy(3,eye(3)*2), ...
    'sphere_geodesic_holonomy:Unit');
verifyError(testCase,@()sphere_geodesic_holonomy(4,eye(3)), ...
    'sphere_geodesic_holonomy:Vertices');
V=eye(3);V(:,2)=-V(:,1);
verifyError(testCase,@()sphere_geodesic_holonomy(3,V), ...
    'sphere_geodesic_holonomy:Antipodal');
verifyError(testCase,@()sphere_geodesic_holonomy(3,[], 'wrong'), ...
    'sphere_geodesic_holonomy:Mode');
verifyError(testCase,@()F(zeros(n-1,1),'notransp'), ...
    'um_spherical_holonomy_apply:Dimension');
verifyError(testCase,@()F(zeros(n,1),'wrong'), ...
    'um_spherical_holonomy_apply:Flag');
end

function test_hyperbolic_area_vs_angle_defect(testCase)
sets={ [0;.45;.35i], [.2+.3i;-.3+.1i;.1-.4i], ...
       [0; .35i; .55], [.65;-.4i;.3+.4i] };
for j=1:numel(sets)
    z=sets{j};
    % Hyperbolic distances in the curvature -1 Poincare disk:
    a=hyp_distance(z(2),z(3));
    b=hyp_distance(z(1),z(3));
    c=hyp_distance(z(1),z(2));
    A=acos(max(-1,min(1, ...
        (cosh(b)*cosh(c)-cosh(a))/(sinh(b)*sinh(c)))));
    B=acos(max(-1,min(1, ...
        (cosh(a)*cosh(c)-cosh(b))/(sinh(a)*sinh(c)))));
    C=acos(max(-1,min(1, ...
        (cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))));
    unsigned=pi-A-B-C;
    H=hyperbolic_triangle_holonomy(2,z);
    angle=atan2(H(2,1),H(1,1));
    oriented=imag(conj(z(2)-z(1))*(z(3)-z(1)));
    signed=sign(oriented)*unsigned;
    verifyLessThan(testCase,abs(angle+signed),2e-11);
    verifyLessThan(testCase,norm(H'*H-eye(2),'fro'),3e-14);
    verifyLessThan(testCase,abs(det(H)-1),3e-14);
    S=hyperbolic_triangle_holonomy(2,z,'spectral');
    verifyEqual(testCase,S.finite.kind,'exact');
    for k=-3:4
        verifyLessThan(testCase,abs(S.finite.moments(k)-trace(H^k)/2),3e-11);
    end
    R=hyperbolic_triangle_holonomy(2,z([1 3 2]));
    verifyLessThan(testCase,norm(R-H','fro'),3e-14);
end
end

function test_hyperbolic_padded_spectral_and_operator(testCase)
for n=[2 3 5 13]
    z=[0;.45;.35i];
    H=hyperbolic_triangle_holonomy(n,z);
    verifyLessThan(testCase,norm(H'*H-eye(n),'fro'),5e-12*n);
    verifyLessThan(testCase,norm(H(3:end,3:end)-eye(n-2),'fro'),5e-12*n);
    S=hyperbolic_triangle_holonomy(n,z,'spectral');
    eigenvalues=S.finite.eigenvalues();
    verifyEqual(testCase,numel(eigenvalues),n);
    verifyLessThan(testCase,max(arrayfun(@(w)min(abs(w-eigenvalues)),eig(H))), ...
        3e-11*n);
    F=hyperbolic_triangle_holonomy(n,z,'operator');
    verifyEqual(testCase,F([],'size'),[n n]);
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(testCase,norm(F(X,'notransp')-H*X,'fro'),3e-11*n);
    verifyLessThan(testCase,norm(F(X,'transp')-H'*X,'fro'),3e-11*n);
end
U=hyperbolic_triangle_holonomy(2,[0;.1;.3]);
verifyLessThan(testCase,norm(U-eye(2),'fro'),3e-14);
end

function test_hyperbolic_validation(testCase)
verifyError(testCase,@()hyperbolic_triangle_holonomy(1),'MATLAB:notGreaterEqual');
verifyError(testCase,@()hyperbolic_triangle_holonomy(2,[1;.1;.2]), ...
    'hyperbolic_triangle_holonomy:Vertices');
verifyError(testCase,@()hyperbolic_triangle_holonomy(2,[0;.2]), ...
    'hyperbolic_triangle_holonomy:Vertices');
verifyError(testCase,@()hyperbolic_triangle_holonomy(2,[],'bad'), ...
    'hyperbolic_triangle_holonomy:Mode');
F=hyperbolic_triangle_holonomy(2,[],'operator');
verifyError(testCase,@()F(ones(1,1),'notransp'), ...
    'um_hyperbolic_holonomy_apply:Dimension');
verifyError(testCase,@()F(ones(2,1),'bad'), ...
    'um_hyperbolic_holonomy_apply:Flag');
end

function test_wilson_default_su_and_unitarity(testCase)
for n=[1 2 3 5 8]
    W=lattice_wilson_loop(n);
    L=lattice_wilson_loop(n,[],'links');
    verifyEqual(testCase,size(L,3),4);
    verifyLessThan(testCase,norm(W'*W-eye(n),'fro'),5e-11*n);
    verifyLessThan(testCase,abs(det(W)-1),5e-10*n);
    for k=1:size(L,3)
        verifyLessThan(testCase,norm(L(:,:,k)'*L(:,:,k)-eye(n),'fro'), ...
            5e-11*n);
        verifyLessThan(testCase,abs(det(L(:,:,k))-1),5e-10*n);
    end
    P=eye(n);
    for k=1:size(L,3),P=P*L(:,:,k);end
    verifyLessThan(testCase,norm(W-P,'fro'),2e-11*n);
    S=lattice_wilson_loop(n,[],'spectral');
    verifyEqual(testCase,S.finite.kind,'unknown');
end
L=lattice_wilson_loop(4,[],'links');
verifyGreaterThan(testCase,norm(L(:,:,1)*L(:,:,2)-L(:,:,2)*L(:,:,1),'fro'),1e-5);
end

function test_wilson_gauge_conjugation_covariance(testCase)
n=4;L=lattice_wilson_loop(n,[],'links');
m=size(L,3);G=zeros(n,n,m);LG=zeros(n,n,m);
for j=1:m,G(:,:,j)=haar_unitary(n,200+j);end
for j=1:m
    next=mod(j,m)+1;
    LG(:,:,j)=G(:,:,j)*L(:,:,j)*G(:,:,next)';
end
W=lattice_wilson_loop(n,L);
WG=lattice_wilson_loop(n,LG);
verifyLessThan(testCase,norm(WG-G(:,:,1)*W*G(:,:,1)','fro'), ...
    5e-11*n);
verifyLessThan(testCase,abs(trace(WG)-trace(W)),5e-11*n);
verifyLessThan(testCase,abs(det(WG)-det(W)),5e-11*n);
z=eig(W);w=eig(WG);
verifyLessThan(testCase,max(arrayfun(@(a)min(abs(a-w)),z)),5e-11*n);
% Orientation reversal maps every link to the adjoint in reverse order.
LR=complex(zeros(size(L)));
for j=1:m,LR(:,:,j)=L(:,:,m-j+1)';end
WR=lattice_wilson_loop(n,LR);
verifyLessThan(testCase,norm(WR-W','fro'),5e-11*n);
end

function test_wilson_operator_arbitrary_u_links_and_validation(testCase)
n=6;m=5;
L=complex(zeros(n,n,m));
for j=1:m,L(:,:,j)=haar_unitary(n,100+j);end
W=lattice_wilson_loop(n,L);
verifyLessThan(testCase,norm(W'*W-eye(n),'fro'),4e-11*n);
verifyLessThan(testCase,abs(det(W)-prod(arrayfun(@(j)det(L(:,:,j)),1:m))), ...
    5e-10*n);
F=lattice_wilson_loop(n,L,'operator');
X=randn(n,3)+1i*randn(n,3);
verifyEqual(testCase,F([],'size'),[n n]);
verifyLessThan(testCase,norm(F(X,'notransp')-W*X,'fro'),3e-11*n);
verifyLessThan(testCase,norm(F(X,'ctransp')-W'*X,'fro'),3e-11*n);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),3e-11*n);
L(:,:,1)=zeros(n);
verifyError(testCase,@()lattice_wilson_loop(n,L),'lattice_wilson_loop:Unitary');
verifyError(testCase,@()lattice_wilson_loop(n,ones(3,3,2)), ...
    'lattice_wilson_loop:Links');
verifyError(testCase,@()lattice_wilson_loop(n,[],'wrong'), ...
    'lattice_wilson_loop:Mode');
verifyError(testCase,@()F(zeros(n-1,1),'notransp'), ...
    'um_wilson_apply:Dimension');
verifyError(testCase,@()F(zeros(n,1),'wrong'), ...
    'um_wilson_apply:Flag');
end

function d=hyp_distance(z,w)
d=2*asinh(abs(z-w)/sqrt((1-abs(z)^2)*(1-abs(w)^2)));
end
