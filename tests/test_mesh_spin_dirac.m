function tests=test_mesh_spin_dirac
%TEST_MESH_SPIN_DIRAC Ye 2018 face-quaternion intrinsic Dirac + area metric.
tests=functiontests(localfunctions);
end

function test_exact_two_triangle_quaternion_spectrum(tc)
M=mesh_spin_dirac([],[],.73,'model');
D=M.D_integrated;H=M.H;N=M.dimension;
verifyEqual(tc,N,8);
verifyTrue(tc,issparse(D));
verifyTrue(tc,issparse(H));
verifyEqual(tc,M.boundary_edge_count,4);
verifyLessThan(tc,norm(full(D-D'),'fro'),3e-14);
verifyLessThan(tc,norm(full(H-H'),'fro'),3e-14);
verifyLessThan(tc,norm(full(D'*M.metric-M.metric*D),'fro'), ...
    1e-13); % D0 is NOT the physical generalized operator when M nonconstant
verifyLessThan(tc,norm(M.U-expm(-.73i*full(H)),'fro'), ...
    1e-13);
verifyLessThan(tc,norm(M.U'*M.U-eye(N),'fro'),2e-13);
verifyLessThan(tc,norm(full(M.D_physical'*M.metric- ...
    M.metric*M.D_physical),'fro'),3e-13);
verifyLessThan(tc,norm(M.V_physical'*M.metric* ...
    M.V_physical-M.metric,'fro'),3e-13);
S=mesh_spin_dirac([],[],.73,'spectral');
verifyEqual(tc,S.finite.kind,'exact');
lam=S.finite.eigenvalues();
verifyEqual(tc,numel(lam),8);
expect=[exp(-.73i*sqrt(2))*ones(4,1); ...
    exp(.73i*sqrt(2))*ones(4,1)];
verifyLessThan(tc,norm(lam-expect),1e-13);
verifyLessThan(tc,max(abs(sort(angle(eig(M.U))) ...
    -sort(angle(lam)))),3e-13);
for k=-4:4
    verifyLessThan(tc,abs(trace(M.U^k)/N ...
        -S.finite.moments(k)),6e-13);
end
end

function test_source_hyperedge_and_nonuniform_face_area(tc)
% A planar two-face mesh has theta=0: each shared hyperedge is
% purely imaginary and has norm equal to physical edge length.
V=[0 0 0;2 0 0;1 1 0;0 1 0];
F=[1 2 3;1 3 4];
M=mesh_spin_dirac(V,F,.82,'model');
verifyEqual(tc,M.boundary_edge_count,4);
verifyLessThan(tc,max(abs(M.edge_signed_angles)),2e-14);
verifyLessThan(tc,norm(full(M.D_integrated(1:4,5:8)) ...
    -full(M.D_integrated(5:8,1:4))','fro'),3e-14);
verifyLessThan(tc,norm(full(M.D_integrated(1:4,5:8))'* ...
    full(M.D_integrated(1:4,5:8)) ...
    -(M.edge_lengths(1)^2/4)*eye(4),'fro'),2e-13);
verifyLessThan(tc,norm(full(M.D_physical'*M.metric ...
    -M.metric*M.D_physical),'fro'),2e-12);
verifyLessThan(tc,norm(full(M.H-M.H'),'fro'),2e-12);
verifyLessThan(tc,norm(M.U'*M.U-eye(8),'fro'),2e-12);
verifyLessThan(tc,norm(M.V_physical'*M.metric ...
    *M.V_physical-M.metric,'fro'),2e-12);
verifyGreaterThan(tc,norm(M.V_physical'*M.V_physical-eye(8),'fro'),.05);
S=mesh_spin_dirac(V,F,.82,'spectral');
verifyEqual(tc,S.finite.kind,'exact');
verifyLessThan(tc,max(arrayfun(@(z)min(abs(eig(M.U)-z)), ...
    S.finite.eigenvalues())),3e-12);
end

function test_curved_tetrahedron_and_local_gauge_covariance(tc)
[V,F]=tetrahedron();
t=.41;
A=mesh_spin_dirac(V,F,t,'model',[],'closed');
verifyEqual(tc,A.boundary_edge_count,0);
verifyEqual(tc,A.dimension,16);
verifyTrue(tc,any(abs(A.edge_signed_angles)>.1));
verifyLessThan(tc,norm(A.U'*A.U-eye(16),'fro'),3e-12);
verifyLessThan(tc,norm(full(A.D_integrated-A.D_integrated'),'fro'), ...
    3e-13);
G=[cos(.2) sin(.2) 0 0; ...
   cos(.3) 0 sin(.3) 0; ...
   cos(.4) 0 0 sin(.4); ...
   1 0 0 0];
B=mesh_spin_dirac(V,F,t,'model',G,'closed');
W=B.gauge_matrix;
verifyLessThan(tc,norm(full(W'*W-speye(16)),'fro'),3e-13);
verifyLessThan(tc,norm(full(B.D_integrated-W'*A.D_integrated*W), ...
    'fro'),3e-13);
verifyLessThan(tc,norm(full(B.H-W'*A.H*W),'fro'),3e-13);
verifyLessThan(tc,norm(B.U-full(W'*A.U*W),'fro'),3e-12);
verifyLessThan(tc,norm(full(B.D_physical'*B.metric ...
    -B.metric*B.D_physical),'fro'),3e-12);
verifyLessThan(tc,norm(B.U'*B.U-eye(16),'fro'),3e-12);
S=mesh_spin_dirac(V,F,t,'spectral',G,'closed');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEqual(tc,S.dimension,16);
end

function test_cayley_rational_approximate_exponential_and_adjoints(tc)
t=.87;
U=mesh_spin_dirac([],[],t);
Y=reshape(sin(.31*(1:24)),8,3) ...
    +1i*reshape(cos(.17*(1:24)),8,3);
errors=zeros(3,1);
for k=1:3
    steps=8*4^(k-1);
    F=mesh_spin_dirac([],[],t,'cayley_operator',[],[],steps);
    verifyEqual(tc,F([],'size'),[8 8]);
    errors(k)=norm(F(eye(8),'notransp')-U,'fro');
    verifyLessThan(tc,norm(F(Y,'notransp')-U*Y,'fro'), ...
        2*errors(k)*norm(Y,'fro'));
    verifyLessThan(tc,norm(F(Y,'ctransp')-U'*Y,'fro'), ...
        2*errors(k)*norm(Y,'fro'));
    verifyLessThan(tc,norm(F(Y,'transp')-U.'*Y,'fro'), ...
        2*errors(k)*norm(Y,'fro'));
    verifyLessThan(tc,norm(F(F(Y,'notransp'),'ctransp')-Y, ...
        'fro'),1e-11*norm(Y,'fro'));
end
verifyGreaterThan(tc,errors(1),errors(2)*10);
verifyGreaterThan(tc,errors(2),errors(3)*10);
verifyLessThan(tc,errors(3),5e-5);
end

function test_orientation_boundary_and_input_rejections(tc)
[V,F]=tetrahedron();
verifyError(tc,@()mesh_spin_dirac([],[],.4,'generator',[],'closed'), ...
    'mesh_spin_dirac:Boundary');
S=mesh_spin_dirac([],[],0,'spectral');
verifyLessThan(tc,norm(S.finite.eigenvalues()-ones(8,1)),1e-14);
Fbad=F;
Fbad(2,[1 2])=Fbad(2,[2 1]);
verifyError(tc,@()mesh_spin_dirac(V,Fbad), ...
    'mesh_spin_dirac:Orientation');
verifyError(tc,@()mesh_spin_dirac([],[],.4,'bad'), ...
    'mesh_spin_dirac:Mode');
verifyError(tc,@()mesh_spin_dirac([],[],.4,'model',ones(2,4)), ...
    'mesh_spin_dirac:Gauges');
verifyError(tc,@()mesh_spin_dirac([],[],.4,'model',[],[],-1), ...
    'mesh_spin_dirac:Steps');
verifyError(tc,@()mesh_spin_dirac([],[],.4,'model',[],'bad'), ...
    'mesh_spin_dirac:Boundary');
Q=mesh_spin_dirac([],[],.4,'cayley_operator');
verifyError(tc,@()Q(ones(7,1),'notransp'), ...
    'um_mesh_spin_dirac_cayley:Dimension');
verifyError(tc,@()Q(ones(8,1),'bad'), ...
    'um_mesh_spin_dirac_cayley:Flag');
end

function [V,F]=tetrahedron()
V=[1 1 1;-1 -1 1;-1 1 -1;1 -1 -1];
F=[1 3 2;1 2 4;2 3 4;3 1 4];
end
