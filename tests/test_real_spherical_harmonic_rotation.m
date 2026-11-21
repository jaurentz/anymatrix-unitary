function tests=test_real_spherical_harmonic_rotation
%TEST_REAL_SPHERICAL_HARMONIC_ROTATION Independent SO(3) real-basis tests.
tests=functiontests(localfunctions);
end

function test_condon_shortley_change_of_basis(tc)
for L=[1 2 3 5 9]
    bases=real_spherical_harmonic_rotation(L,.3,.7,.4,'basis');
    A=real_spherical_harmonic_rotation(L,.3,.7,.4);
    B=spherical_harmonic_rotation(L,.3,.7,.4);
    verifyTrue(tc,isreal(A));
    verifyTrue(tc,issparse(A));
    verifyEqual(tc,size(A),[L^2 L^2]);
    verifyLessThan(tc,norm(full(A'*A-speye(L^2)),'fro'), ...
        8e-11*L^2);
    verifyLessThan(tc,norm(full(B'*B-speye(L^2)),'fro'), ...
        8e-11*L^2);
    offset=0;
    for ell=0:L-1
        inds=offset+(1:2*ell+1);
        C=bases{ell+1};
        Dc=wigner_rotation(ell,.3,.7,.4);
        Dr=full(A(inds,inds));
        verifyLessThan(tc,norm(C'*C-eye(2*ell+1),'fro'), ...
            4e-14*(2*ell+1));
        verifyLessThan(tc,norm(Dr-C'*Dc*C,'fro'), ...
            6e-12*(2*ell+1));
        verifyLessThan(tc,norm(imag(C'*Dc*C),'fro'), ...
            2e-12*(2*ell+1));
        % The complex convention is not the entrywise real part.
        if ell>0
            verifyGreaterThan(tc,norm(Dr-real(Dc),'fro'),1e-3);
        end
        offset=offset+2*ell+1;
    end
end
end

function test_real_so3_group_law_and_singular_euler_limits(tc)
a1=.3;b1=.7;g1=-.4;
a2=-.6;b2=1.1;g2=.2;
R=physical_zyz(a1,b1,g1)*physical_zyz(a2,b2,g2);
bc=acos(max(-1,min(1,R(3,3))));
ac=atan2(R(2,3),R(1,3));
gc=atan2(R(3,2),-R(3,1));
for L=[1 2 3 5 7]
    Q1=real_spherical_harmonic_rotation(L,a1,b1,g1);
    Q2=real_spherical_harmonic_rotation(L,a2,b2,g2);
    Q3=real_spherical_harmonic_rotation(L,ac,bc,gc);
    verifyLessThan(tc,norm(full(Q1*Q2-Q3),'fro'), ...
        7e-10*L^2);
    verifyLessThan(tc,norm(full(Q1'*Q1-speye(L^2)),'fro'), ...
        8e-11*L^2);
    for beta=[0 pi -pi]
        O=real_spherical_harmonic_rotation(L,.4,beta,.9);
        verifyTrue(tc,isreal(O));
        verifyLessThan(tc,norm(full(O'*O-speye(L^2)),'fro'), ...
            8e-11*L^2);
    end
end
end

function test_exact_multiplicity_and_trace_moments(tc)
for L=[1 2 3 6 11]
    a=.3;b=.7;g=.4;
    U=real_spherical_harmonic_rotation(L,a,b,g);
    S=real_spherical_harmonic_rotation(L,a,b,g,'spectral');
    verifyEqual(tc,S.family,'real_spherical_harmonic_rotation');
    verifyEqual(tc,S.finite.kind,'exact');
    z=S.finite.eigenvalues();
    verifyEqual(tc,numel(z),L^2);
    delta=2*acos(cos(b/2)*cos((a+g)/2));
    target=[];
    for m=-(L-1):L-1
        target=[target;repmat(exp(-1i*m*delta),L-abs(m),1)]; %#ok<AGROW>
    end
    verifyEqual(tc,numel(target),L^2);
    for j=1:numel(z)
        verifyLessThan(tc,min(abs(target-z(j))),1e-10*L);
    end
    for k=-3:3
        observed=trace(full(U)^k)/L^2;
        expect=sum(target.^k)/L^2;
        verifyLessThan(tc,abs(observed-expect),5e-10*L^2);
        verifyLessThan(tc,abs(S.finite.moments(k)-expect),5e-12*L^2);
    end
end
end

function test_complex_and_real_blockwise_operators(tc)
for L=[1 2 5 9 16]
    A=spherical_harmonic_rotation(L,.19,.67,-.34);
    B=real_spherical_harmonic_rotation(L,.19,.67,-.34);
    FA=spherical_harmonic_rotation(L,.19,.67,-.34,'operator');
    FB=real_spherical_harmonic_rotation(L,.19,.67,-.34,'operator');
    N=L^2;
    verifyEqual(tc,FA([],'size'),[N N]);
    verifyEqual(tc,FB([],'size'),[N N]);
    X=randn(N,3)+1i*randn(N,3);
    ops={FA,FB};matrices={A,B};
    for caseIndex=1:2
        op=ops{caseIndex};matrix=matrices{caseIndex};
        verifyLessThan(tc,norm(op(X,'notransp')-matrix*X,'fro'), ...
            2e-10*N);
        verifyLessThan(tc,norm(op(X,'transp')-matrix.'*X,'fro'), ...
            2e-10*N);
        verifyLessThan(tc,norm(op(X,'ctransp')-matrix'*X,'fro'), ...
            2e-10*N);
        verifyLessThan(tc,norm(op(op(X,'notransp'),'ctransp')-X,'fro'), ...
            2e-10*N);
    end
end
end

function test_bounded_explicit_and_invalid_modes(tc)
verifyError(tc,@()real_spherical_harmonic_rotation(73), ...
    'real_spherical_harmonic_rotation:ExplicitSize');
verifyError(tc,@()real_spherical_harmonic_rotation(3,.1,.2,.3,'bad'), ...
    'real_spherical_harmonic_rotation:Mode');
verifyError(tc,@()spherical_harmonic_rotation(3,.1,.2,.3,'bad'), ...
    'spherical_harmonic_rotation:Mode');
F=real_spherical_harmonic_rotation(6,.1,.2,.3,'operator');
verifyError(tc,@()F(zeros(35,1),'notransp'), ...
    'um_wigner_direct_sum_apply:Dimension');
verifyError(tc,@()F(zeros(36,1),'bad'), ...
    'um_wigner_direct_sum_apply:Flag');
end

function R=physical_zyz(a,b,g)
Za=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1];
Yb=[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
Zg=[cos(g) -sin(g) 0;sin(g) cos(g) 0;0 0 1];
R=Za*Yb*Zg;
end
