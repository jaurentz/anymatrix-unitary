function tests=test_polar_galerkin_stages
%TEST_POLAR_GALERKIN_STAGES Source equations (24b,c) only; NOT full VPD.
tests=functiontests(localfunctions);
end

function test_source_table1_gauss_legendre_collocation_coefficients(tc)
for s=1:3
    [A,b,c,order]=um_polar_galerkin_coefficients(s);
    verifyEqual(tc,size(A),[s s]);
    verifyEqual(tc,numel(b),s);
    verifyEqual(tc,numel(c),s);
    verifyEqual(tc,order,2*s);
    verifyLessThan(tc,norm(A*ones(s,1)-c),2e-14);
    verifyLessThan(tc,abs(sum(b)-1),2e-14);
    verifyLessThan(tc,norm(diag(b)*A+A'*diag(b)-b*b','fro'), ...
        3e-14);
end
verifyError(tc,@()um_polar_galerkin_coefficients(4), ...
    'um_polar_galerkin_coefficients:Stages');
end

function test_fixed_point_equations_and_SO3_SO4(tc)
for n=[3 4]
    W=reshape(sin(.13*(1:n*n)),n,n);
    Omega=.16*(W-W');
    G0=expm(.3*Omega);
    for s=1:3
        [G1,U,info]=um_polar_galerkin_stages( ...
            G0,Omega,.17,s,1e-12,250);
        verifyTrue(tc,info.iterations>0);
        verifyFalse(tc,info.is_full_variational_integrator);
        verifyFalse(tc,info.has_momentum_or_symplectic_update);
        verifyLessThan(tc,info.stage_residual,1e-12);
        verifyLessThan(tc,info.maximum_orthogonality_defect, ...
            2e-12*n);
        verifyLessThan(tc,norm(G1'*G1-eye(n),'fro'),2e-12*n);
        verifyLessThan(tc,abs(det(G1)-1),3e-12);
        for i=1:s
            verifyLessThan(tc,norm(U(:,:,i)'*U(:,:,i)-eye(n),'fro'), ...
                2e-12*n);
            A=G0;
            for j=1:s
                A=A+.17*info.A(i,j)*U(:,:,j)*Omega;
            end
            Q=independent_polar(A);
            verifyLessThan(tc,norm(U(:,:,i)-Q,'fro'),5e-12*n);
        end
        A=G0;
        for j=1:s
            A=A+.17*info.b(j)*U(:,:,j)*Omega;
        end
        verifyLessThan(tc,norm(G1-independent_polar(A),'fro'), ...
            3e-12*n);
    end
end
end

function test_left_equivariance_noncommuting_given_stage_omegas(tc)
n=4;
W=reshape(cos(.19*(1:n*n)),n,n);
K=.22*(W-W');
G0=expm(.23*K);
S=expm(.31*K);
Om=zeros(n,n,3);
for j=1:3
    B=reshape(sin(.11*(1:n*n)+.4*j),n,n);
    Om(:,:,j)=.15*(B-B');
end
for s=1:3
    [Q,Ui]=um_polar_galerkin_stages(G0,Om(:,:,1:s),.08,s);
    [R,Vi]=um_polar_galerkin_stages(S*G0,Om(:,:,1:s),.08,s);
    verifyLessThan(tc,norm(R-S*Q,'fro'),5e-11*n);
    for j=1:s
        verifyLessThan(tc,norm(Vi(:,:,j)-S*Ui(:,:,j),'fro'), ...
            5e-11*n);
    end
end
end

function test_commuting_omega_rotation_local_order(tc)
% Prescribed CONSTANT Omega is an exactly solvable reference flow.
% This validates the stage equations' numerical consistency;
% it is NOT a symplectic-order test for the missing momentum update.
W=[0 -1.2 .3;1.2 0 -.4;-.3 .4 0];
for s=1:3
    E=zeros(3,1);
    for k=1:3
        h=.2/(2^(k-1));
        [G1,~,info]=um_polar_galerkin_stages(eye(3),W,h,s);
        E(k)=norm(G1-expm(h*W),'fro');
        verifyLessThan(tc,info.stage_residual,5e-12);
    end
    ratios=[5 20 55];
    verifyGreaterThan(tc,E(1),E(2)*ratios(s));
    verifyGreaterThan(tc,E(2),E(3)*ratios(s));
end
end

function test_nonlinear_solver_and_domain_failures(tc)
W=[0 -1 0;1 0 0;0 0 0];
verifyError(tc,@()um_polar_galerkin_stages(eye(3),W,.2,2, ...
    1e-12,1),'um_polar_galerkin_stages:NoConvergence');
verifyError(tc,@()um_polar_galerkin_stages(2*eye(3),W,.2), ...
    'um_polar_galerkin_stages:Configuration');
verifyError(tc,@()um_polar_galerkin_stages(diag([-1 1 1]),W,.2), ...
    'um_polar_galerkin_stages:Configuration');
verifyError(tc,@()um_polar_galerkin_stages(eye(3),eye(3),.2), ...
    'um_polar_galerkin_stages:Omega');
verifyError(tc,@()um_polar_galerkin_stages(eye(3),W,Inf), ...
    'um_polar_galerkin_stages:Step');
end

function Q=independent_polar(A)
[U,~,V]=svd(A);
Q=U*V';
end
