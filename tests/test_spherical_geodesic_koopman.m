function tests=test_spherical_geodesic_koopman
%TEST_SPHERICAL_GEODESIC_KOOPMAN Exact Koopman-Wigner equivalence (#61).
tests=functiontests(localfunctions);
end

function test_round_sphere_frame_geodesic_and_degree_one(tc)
% The unit tangent frame g=[x,v,x cross v] follows right translation.
B=[1/sqrt(2),0,-1/sqrt(2); ...
   -1i/sqrt(2),0,-1i/sqrt(2); ...
   0,1,0];
verifyLessThan(tc,norm(B'*B-eye(3),'fro'),1e-14);
for t=[0 .37 -1.2 pi]
    R=[cos(t),-sin(t),0;sin(t),cos(t),0;0,0,1];
    G0=[cos(.28),-sin(.28),0;sin(.28),cos(.28),0;0,0,1] ...
       *[cos(.35),0,sin(.35);0,1,0;-sin(.35),0,cos(.35)];
    x=G0(:,1);v=G0(:,2);
    xt=cos(t)*x+sin(t)*v;
    vt=-sin(t)*x+cos(t)*v;
    Gt=[xt,vt,cross(xt,vt)];
    verifyLessThan(tc,norm(Gt-G0*R,'fro'),2e-14);
    verifyLessThan(tc,norm(B'*R*B-wigner_rotation(1,t,0,0), ...
        'fro'),2e-13);
end
end

function test_exact_finite_wigner_koopman_blocks(tc)
for ell=0:6
    d=2*ell+1;
    m=(-ell:ell)';
    for t=[0 .24 -1.13 2*pi/7]
        % Fixed-left-index Peter-Weyl subspace: right translation acts
        % by the existing degree-ell Wigner representation.
        W=wigner_rotation(ell,t,0,0);
        expected=diag(exp(-1i*m*t));
        verifyEqual(tc,size(W),[d d]);
        verifyLessThan(tc,norm(W-expected,'fro'),2e-12*d);
        verifyLessThan(tc,norm(W'*W-eye(d),'fro'),2e-12*d);
        % Full degree-ell matrix coefficient space repeats this
        % block d times, giving each phase multiplicity d.
        K=kron(speye(d),sparse(W));
        verifyLessThan(tc,norm(full(K'*K-speye(d*d)),'fro'), ...
            3e-12*d);
        lam=repmat(exp(-1i*m*t),d,1);
        for k=0:4
            verifyLessThan(tc,abs(trace(full(K)^k) ...
                -sum(lam.^k)),3e-11*d*d);
        end
        verifyLessThan(tc,norm(full(K-kron(speye(d), ...
            sparse(expected))),'fro'),2e-12*d*d);
    end
end
end
