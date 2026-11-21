function tests=test_moser_veselov
%TEST_MOSER_VESELOV Independent discrete rigid-body equation benchmarks.
tests=functiontests(localfunctions);
end

function test_published_SO3_planar_equation_specialization(tc)
% Moser--Veselov / Cardoso--Miraldo equation F'*J-JF=h*M,
% specialized to a principal-axis rotation in SO(3), with an
% independent exact solution of its only nonzero 2x2 skew block.
h=.24;phi=.18;
J=diag([.8 1.3 2]);
R=[cos(phi) -sin(phi);sin(phi) cos(phi)];
Ftrue=blkdiag(R,1);
M=zeros(3);
M(1,2)=(J(1,1)+J(2,2))*sin(phi)/h;
M(2,1)=-M(1,2);
F=moser_veselov(3,h,J,M);
verifyLessThan(tc,norm(F-Ftrue,'fro'),2e-11);
verifyLessThan(tc,norm(F'*J-J*F-h*M,'fro'),3e-12);
verifyLessThan(tc,norm(F'*F-eye(3),'fro'),2e-12);
verifyEqual(tc,det(F),1,'AbsTol',2e-12);
verifyLessThan(tc,min(abs(eig(F)-1)),2e-12);
end

function test_known_SO4_SO6_manufactured_rotations(tc)
% Independent rotations are prescribed ONLY to manufacture
% admissible test momentum; generator receives J and M (not F).
for n=[4 5 6]
    h=.16;d=(1:n)';
    J=diag(1+.2*d)+.035*(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1));
    A=zeros(n);
    for j=1:n-1,A(j,j+1)=.12*sin(.37*j);end
    for j=1:n-2,A(j,j+2)=.08*cos(.29*j);end
    A=A-A';
    Ftrue=expm(-h*A);
    M=(Ftrue'*J-J*Ftrue)/h;
    F=moser_veselov(n,h,J,M);
    verifyLessThan(tc,norm(F-Ftrue,'fro'),3e-11*n);
    verifyLessThan(tc,norm(F'*J-J*F-h*M,'fro'),3e-11*n);
    verifyLessThan(tc,norm(F'*F-eye(n),'fro'),3e-12*n);
    verifyEqual(tc,det(F),1,'AbsTol',3e-12*n);
end
end

function test_full_SO_n_nontrivial_and_inertia(tc)
for n=[3 4 5 8 12]
    F=moser_veselov(n);
    verifyEqual(tc,size(F),[n n]);
    verifyLessThan(tc,norm(F'*F-eye(n),'fro'),2e-11*n);
    verifyEqual(tc,det(F),1,'AbsTol',2e-11*n);
    verifyGreaterThan(tc,norm(F-eye(n),'fro'),1e-4);
    D=moser_veselov(n,[],[],[],1,'model');
    verifyLessThan(tc,norm(D.F-F,'fro'),2e-12*n);
    verifyLessThan(tc,D.residual,2e-10);
    verifyLessThan(tc,D.spatial_momentum_residual,2e-10);
    verifyLessThan(tc,norm(D.M+D.M','fro'),1e-13);
    verifyGreaterThan(tc,min(eig(D.J)),0);
    verifyLessThan(tc,norm(D.M_next-F*D.M*F','fro'),3e-12);
    verifyLessThan(tc,norm(D.Q_next-F','fro'),3e-12);
    vals=eig(F);
    verifyLessThan(tc,max(abs(abs(vals)-1)),3e-11);
    verifyLessThan(tc,max(arrayfun(@(z)min(abs(vals-conj(z))),vals)), ...
        1e-9);
    if mod(n,2)==1
        verifyLessThan(tc,min(abs(vals-1)),1e-9);
    end
end
end

function test_discrete_Euler_Arnold_orbit_and_spatial_momentum(tc)
n=6;steps=12;
B=moser_veselov(n,.14,[],[],steps,'trajectory');
verifyEqual(tc,size(B.F),[n n steps]);
verifyEqual(tc,size(B.Q),[n n steps+1]);
verifyEqual(tc,size(B.M),[n n steps+1]);
verifyLessThan(tc,norm(B.Q(:,:,1)-eye(n),'fro'),1e-14);
verifyLessThan(tc,max(B.residuals),1e-10);
M0=B.M(:,:,1);
fro0=norm(M0,'fro');
for j=1:steps
    F=B.F(:,:,j);Q=B.Q(:,:,j);Qnext=B.Q(:,:,j+1);
    M=B.M(:,:,j);Mnext=B.M(:,:,j+1);
    verifyLessThan(tc,norm(F'*F-eye(n),'fro'),5e-11*n);
    verifyEqual(tc,det(F),1,'AbsTol',5e-11);
    verifyLessThan(tc,norm(F'*B.J-B.J*F-B.h*M,'fro'), ...
        3e-11*n);
    verifyLessThan(tc,norm(Mnext-F*M*F','fro'),2e-11*n);
    verifyLessThan(tc,norm(Qnext-Q*F','fro'),2e-11*n);
    verifyLessThan(tc,norm(Qnext'*Qnext-eye(n),'fro'),2e-11*n);
    verifyLessThan(tc,norm(Qnext*Mnext*Qnext'-M0,'fro'), ...
        3e-11*n);
    verifyEqual(tc,norm(Mnext,'fro'),fro0,'AbsTol',2e-11);
    verifyLessThan(tc,norm(Mnext+Mnext','fro'),2e-11);
end
end

function test_zero_momentum_identity_and_finite_spectral_reference(tc)
for n=[3 4 5 8]
    F=moser_veselov(n,.45,[],zeros(n));
    verifyLessThan(tc,norm(F-eye(n),'fro'),1e-13);
    S=moser_veselov(n,.45,[],zeros(n),1,'spectral');
    verifyEqual(tc,S.dimension,n);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.finite.eigenvalues(),ones(n,1));
    verifyEqual(tc,S.finite.atoms.angles,0);
    verifyEqual(tc,S.finite.atoms.weights,1);
    verifyEqual(tc,S.finite.moments([-3 0 1 2]),ones(1,4));
end
S=moser_veselov(4,.12,[],[],1,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
end

function test_validations_no_nonorthogonal_fallback(tc)
[A,properties]=moser_veselov();
verifyEmpty(tc,A);
for prop={'orthogonal','unitary','real','square'}
    verifyTrue(tc,any(strcmp(prop{1},properties)));
end
verifyError(tc,@()moser_veselov(2), ...
    'MATLAB:expectedGreaterThanOrEqual');
verifyError(tc,@()moser_veselov(4,-1), ...
    'MATLAB:expectedPositive');
verifyError(tc,@()moser_veselov(4,.15,diag([1 1 -2 2])), ...
    'moser_veselov:Inertia');
verifyError(tc,@()moser_veselov(4,.15, ...
    [1 1 0 0;0 1 0 0;0 0 1 0;0 0 0 1]), ...
    'moser_veselov:Inertia');
verifyError(tc,@()moser_veselov(4,.15,eye(4),ones(4)), ...
    'moser_veselov:Momentum');
verifyError(tc,@()moser_veselov(4,.15,eye(4), ...
    zeros(3,3)),'moser_veselov:Momentum');
verifyError(tc,@()moser_veselov(4,.15,eye(4),zeros(4),1, ...
    'nonsense'),'moser_veselov:Mode');
% A necessary norm bound for F'*J-JF=h*M rejects impossible data.
M=zeros(4);M(1,2)=5;M(2,1)=-5;
verifyError(tc,@()moser_veselov(4,1,eye(4),M), ...
    'um_moser_veselov_solve:NoSolution');
end
