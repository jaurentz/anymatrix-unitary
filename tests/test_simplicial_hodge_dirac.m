function tests=test_simplicial_hodge_dirac
%TEST_SIMPLICIAL_HODGE_DIRAC Metric/chiral/topological finite unitary.
tests=functiontests(localfunctions);
end

function test_closed_sphere_metric_chiral_and_betti(tc)
M=simplicial_hodge_dirac([],[],.67,'model');
N=M.dimension;
verifyEqual(tc,M.cochain_dimensions,[6 12 8]);
verifyEqual(tc,N,26);
verifyEqual(tc,M.betti,[1 0 1]);
verifyEqual(tc,M.euler_characteristic,2);
verifyEqual(tc,M.boundary_edge_count,0);
verifyTrue(tc,issparse(M.H));
verifyLessThan(tc,norm(full(M.D1*M.D0),'fro'),1e-14);
verifyLessThan(tc,norm(full(M.d*M.d),'fro'),1e-14);
verifyLessThan(tc,norm(full(M.H-M.H'),'fro'),1e-13);
verifyLessThan(tc,norm(M.U'*M.U-eye(N),'fro'),2e-11*N);
verifyLessThan(tc,norm(M.U-M.U.','fro'),2e-11*N);
verifyLessThan(tc,norm(M.U-expm(-1i*.67*full(M.H)),'fro'), ...
    2e-11*N);
verifyLessThan(tc,norm(full(M.D_physical'*M.metric- ...
    M.metric*M.D_physical),'fro'),3e-12*N);
verifyLessThan(tc,norm(full(M.V_physical'*M.metric* ...
    M.V_physical-M.metric),'fro'),4e-11*N);
verifyGreaterThan(tc,norm(full(M.V_physical'* ...
    M.V_physical-eye(N)),'fro'),.01);
G=diag(M.chiral_grading);
verifyLessThan(tc,norm(G*full(M.H)*G+full(M.H),'fro'), ...
    2e-12*N);
verifyLessThan(tc,norm(G*M.U*G-M.U','fro'),2e-11*N);
ev=eig(full(M.H));
verifyEqual(tc,sum(abs(ev)<2e-10),2);
verifyLessThan(tc,norm(sort(ev)+flipud(sort(ev))),3e-10*N);
S=simplicial_hodge_dirac([],[],.67,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEqual(tc,S.guaranteed_unit_eigenvalues,2);
verifyEqual(tc,S.known_zero_modes,2);
verifyEqual(tc,S.betti,[1 0 1]);
% Topological zero-mode phase count is a LOWER bound only:
% at t=0 every cochain state has phase 1, not merely the Betti modes.
S0=simplicial_hodge_dirac([],[],0,'spectral');
verifyEqual(tc,S0.guaranteed_unit_eigenvalues,2);
U0=simplicial_hodge_dirac([],[],0);
verifyLessThan(tc,norm(U0-eye(N),'fro'),2e-13);
U=simplicial_hodge_dirac([],[],.67);
verifyLessThan(tc,norm(U-M.U,'fro'),2e-13);
end

function test_disk_and_torus_topological_zero_modes(tc)
V=[0 0 0;1 0 0;0 1 0];
F=[1 2 3];
M=simplicial_hodge_dirac(V,F,.43,'generator');
verifyEqual(tc,M.cochain_dimensions,[3 3 1]);
verifyEqual(tc,M.betti,[1 0 0]);
verifyEqual(tc,M.boundary_edge_count,3);
verifyEqual(tc,M.euler_characteristic,1);
verifyEqual(tc,sum(abs(eig(full(M.H)))<2e-10),1);
V=[0 0 0;1 0 0;1 1 0;0 1 0];
F=[1 2 3;1 3 4];
M=simplicial_hodge_dirac(V,F,.43,'generator');
verifyEqual(tc,M.betti,[1 0 0]);
verifyEqual(tc,M.boundary_edge_count,4);
[V,F]=embedded_torus(4,4);
M=simplicial_hodge_dirac(V,F,.43,'generator');
verifyEqual(tc,M.betti,[1 2 1]);
verifyEqual(tc,M.euler_characteristic,0);
verifyEqual(tc,M.boundary_edge_count,0);
verifyLessThan(tc,norm(full(M.D1*M.D0),'fro'),1e-13);
verifyEqual(tc,sum(abs(eig(full(M.H)))<5e-10),4);
S=simplicial_hodge_dirac(V,F,.43,'spectral');
verifyEqual(tc,S.guaranteed_unit_eigenvalues,4);
verifyEqual(tc,S.finite.kind,'unknown');
end

function test_orientation_change_rigid_motion_and_scale(tc)
M=simplicial_hodge_dirac([],[],.37,'generator');
V=M.vertices;F=M.faces;
G=simplicial_hodge_dirac(V,F(:,[1 3 2]),.37,'generator');
d=M.cochain_dimensions;
T=diag([ones(d(1)+d(2),1);-ones(d(3),1)]);
verifyLessThan(tc,norm(full(G.H-T*M.H*T),'fro'),2e-13);
verifyEqual(tc,M.betti,G.betti);
angle=.47;
R=[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
Vrot=(R*V.').';
Rot=simplicial_hodge_dirac(Vrot,F,.37,'generator');
verifyLessThan(tc,norm(full(Rot.H-M.H),'fro'),2e-11);
scale=2.31;
Big=simplicial_hodge_dirac(scale*V,F,.37,'generator');
verifyLessThan(tc,norm(full(Big.H-M.H/scale),'fro'),2e-11);
verifyEqual(tc,Big.betti,M.betti);
end

function test_sparse_cayley_unitary_approximation_and_adjoint(tc)
t=.91;
U=simplicial_hodge_dirac([],[],t);
err=zeros(1,4);
for k=1:4
    steps=2^(k+1);
    C=simplicial_hodge_dirac([],[],t,'cayley_operator',steps);
    F=C.apply;
    n=C.dimension;
    Q=F(eye(n),'notransp');
    verifyTrue(tc,contains(C.approximation,'not_exact_expm'));
    verifyEqual(tc,F([],'size'),[n n]);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),2e-11*n);
    verifyLessThan(tc,norm(Q-Q.','fro'),2e-11*n);
    X=reshape(sin(.11*(1:3*n)),n,3) ...
        +1i*reshape(cos(.23*(1:3*n)),n,3);
    verifyLessThan(tc,norm(F(X,'notransp')-Q*X,'fro'), ...
        2e-11*n);
    verifyLessThan(tc,norm(F(X,'ctransp')-Q'*X,'fro'), ...
        2e-11*n);
    verifyLessThan(tc,norm(F(X,'transp')-Q.'*X,'fro'), ...
        2e-11*n);
    verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X, ...
        'fro'),3e-11*n);
    err(k)=norm(Q-U,'fro');
end
verifyGreaterThan(tc,min(err(1:3)./err(2:4)),3.1);
verifyGreaterThan(tc,err(end),1e-9);
end

function test_sparse_mesh_growth_and_validation(tc)
[empty,props]=simplicial_hodge_dirac();
verifyEmpty(tc,empty);
verifyTrue(tc,any(strcmp(props,'unitary')));
[V,F]=embedded_torus(8,8);
M=simplicial_hodge_dirac(V,F,1,'generator');
verifyEqual(tc,M.dimension,384);
verifyTrue(tc,issparse(M.H));
verifyLessThan(tc,nnz(M.H),40*M.dimension);
verifyError(tc,@()simplicial_hodge_dirac(V,F), ...
    'simplicial_hodge_dirac:DenseOrder');
verifyError(tc,@()simplicial_hodge_dirac(V,F,1,'wrong'), ...
    'simplicial_hodge_dirac:Mode');
verifyError(tc,@()simplicial_hodge_dirac(V,F,1,'generator',1.5), ...
    'simplicial_hodge_dirac:Steps');
C=simplicial_hodge_dirac([],[],1,'cayley_operator');
verifyError(tc,@()C.apply(zeros(25,1),'notransp'), ...
    'um_simplicial_hodge_cayley:Dimension');
verifyError(tc,@()C.apply(zeros(26,1),'wrong'), ...
    'um_simplicial_hodge_cayley:Flag');
V=[0 0 0;1 0 0;0 1 0;0 0 1];
verifyError(tc,@()simplicial_hodge_dirac(V,[1 2 3;1 2 4]), ...
    'simplicial_hodge_dirac:Orientation');
verifyError(tc,@()simplicial_hodge_dirac(V,[1 2 3;1 2 3]), ...
    'simplicial_hodge_dirac:Topology');
verifyError(tc,@()simplicial_hodge_dirac(V,[1 2 3]), ...
    'simplicial_hodge_dirac:Topology');
end

function [V,F]=embedded_torus(nu,nv)
% Independent explicitly triangulated periodic torus; χ=0, genus=1.
V=zeros(nu*nv,3);
for j=0:nv-1
    th=2*pi*j/nv;
    for i=0:nu-1
        ph=2*pi*i/nu;
        idx=i+1+nu*j;
        V(idx,:)=[(3+.65*cos(th))*cos(ph), ...
            (3+.65*cos(th))*sin(ph),.65*sin(th)];
    end
end
F=zeros(2*nu*nv,3);
k=0;
for j=0:nv-1
    jp=mod(j+1,nv);
    for i=0:nu-1
        ip=mod(i+1,nu);
        a=i+1+nu*j;b=ip+1+nu*j;
        c=i+1+nu*jp;d=ip+1+nu*jp;
        k=k+1;F(k,:)=[a b d];
        k=k+1;F(k,:)=[a d c];
    end
end
end
