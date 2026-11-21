function tests=test_unitary_hessenberg_pencil
%TEST_UNITARY_HESSENBERG_PENCIL Independent strict pencil pole checks.
tests=functiontests(localfunctions);
end

function test_strict_core_factors_projective_poles_and_spectra(tc)
for n=[2 3 6 11 20]
    j=(1:n-1)';
    poles=(.27+.19*j).*exp(1i*.27*j);
    for seed=[1 7]
        P=unitary_hessenberg_pencil(n,poles,seed,'complex','pencil');
        U=unitary_hessenberg_pencil(n,poles,seed,'complex');
        D=unitary_hessenberg_pencil(n,poles,seed,'complex','parameters');
        verifyTrue(tc,P.is_proper);
        verifyTrue(tc,all(P.proper));
        verifyEqual(tc,size(P.A),[n n]);
        verifyEqual(tc,size(P.B),[n n]);
        verifyLessThan(tc,norm(P.A'*P.A-eye(n),'fro'),2e-10*n);
        verifyLessThan(tc,norm(P.B'*P.B-eye(n),'fro'),2e-10*n);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),3e-10*n);
        verifyLessThan(tc,norm(P.B'*P.A-U,'fro'),2e-12*n);
        verifyLessThan(tc,norm(tril(P.A,-2),'fro'),2e-12*n);
        verifyLessThan(tc,norm(tril(P.B,-2),'fro'),2e-12*n);
        actual=diag(P.A,-1)./diag(P.B,-1);
        verifyLessThan(tc,norm(actual-poles),5e-9*n);
        verifyLessThan(tc,norm(P.poles-actual),5e-9*n);
        verifyLessThan(tc,norm(D.achieved_poles-poles),5e-9*n);
        verifyEqual(tc,D.core_count,2*(n-1));
        verifyEqual(tc,D.dimension,n);
        verifyLessThan(tc,norm(D.rhoA./D.rhoB-abs(poles)),2e-8*n);
        generalized=eig(P.A,P.B);
        ordinary=eig(U);
        for k=1:n
            verifyLessThan(tc,min(abs(generalized-ordinary(k))), ...
                6e-9*n);
        end
        verifyLessThan(tc,max(abs(abs(generalized)-1)),4e-10*n);
        F=unitary_hessenberg_pencil(n,poles,seed,'complex','operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),3e-10*n);
        verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),3e-10*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),3e-10*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            4e-10*n);
    end
end
end

function test_zero_infinite_and_near_reducible(tc)
poles=[.0001; 0; Inf; 2e3; -.9; -Inf; 1.2+.35i];
n=numel(poles)+1;
P=unitary_hessenberg_pencil(n,poles,23,'complex','pencil');
D=unitary_hessenberg_pencil(n,poles,23,'complex','parameters');
verifyFalse(tc,P.is_proper);
verifyEqual(tc,P.proper,[true;false;false;true;true;false;true]);
verifyEqual(tc,nnz(diag(P.A,-1)==0),1);
verifyEqual(tc,nnz(diag(P.B,-1)==0),2);
verifyEqual(tc,P.poles(2),0);
verifyTrue(tc,isinf(P.poles(3)));
verifyTrue(tc,isinf(P.poles(6)));
idx=[1 4 5 7];
observed=diag(P.A,-1)./diag(P.B,-1);
verifyLessThan(tc,norm(observed(idx)-poles(idx)),1e-5);
verifyTrue(tc,all(isfinite(D.achieved_poles(idx))));
verifyLessThan(tc,norm(P.A'*P.A-eye(n),'fro'),2e-10*n);
verifyLessThan(tc,norm(P.B'*P.B-eye(n),'fro'),2e-10*n);
verifyLessThan(tc,norm(P.U'*P.U-eye(n),'fro'),2e-10*n);
end

function test_real_orthogonal_and_conjugate_spectrum(tc)
n=13;
poles=[-2.3; .3; -1.2; .75; -5; 2.1; ...
    .6; -1.7; Inf; 0; 4.2; -3.4];
P=unitary_hessenberg_pencil(n,poles,42,'real','pencil');
verifyTrue(tc,isreal(P.A));
verifyTrue(tc,isreal(P.B));
verifyTrue(tc,isreal(P.U));
verifyLessThan(tc,norm(P.A'*P.A-eye(n),'fro'),2e-10*n);
verifyLessThan(tc,norm(P.B'*P.B-eye(n),'fro'),2e-10*n);
verifyLessThan(tc,norm(P.U'*P.U-eye(n),'fro'),2e-10*n);
lam=eig(P.U);
for k=1:n
    verifyLessThan(tc,min(abs(lam-conj(lam(k)))),3e-10*n);
end
Q=unitary_hessenberg_pencil(n,poles,42,'real','pencil');
verifyEqual(tc,P.A,Q.A,'AbsTol',0);
verifyEqual(tc,P.B,Q.B,'AbsTol',0);
end

function test_large_operator_and_seed_state(tc)
n=4096;
pole=(.4+.8*sin((1:n-1)'/23).^2);
F=unitary_hessenberg_pencil(n,pole,8,'real','operator');
X=randn(n,2)+1i*randn(n,2);
verifyEqual(tc,F([],'size'),[n n]);
verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'),4e-10);
S=unitary_hessenberg_pencil(7,[],4,[],'spectral');
verifyEqual(tc,S.family,'unitary_hessenberg_pencil');
verifyEqual(tc,S.finite.kind,'unknown');
before=rng;
unitary_hessenberg_pencil(20,[],1,'complex','parameters');
verifyEqual(tc,rng,before);
end

function test_reject_unphysical_parameter_claims(tc)
verifyError(tc,@()unitary_hessenberg_pencil(5,[.5;1;2;1e-5]), ...
    'unitary_hessenberg_pencil:ResolvedPoles');
verifyError(tc,@()unitary_hessenberg_pencil(4,[1;1i;2],1,'real'), ...
    'unitary_hessenberg_pencil:Real');
verifyError(tc,@()unitary_hessenberg_pencil(4,[1;NaN;2]), ...
    'unitary_hessenberg_pencil:Poles');
verifyError(tc,@()unitary_hessenberg_pencil(3,[],1,'complex','bad'), ...
    'unitary_hessenberg_pencil:Mode');
verifyError(tc,@()unitary_hessenberg_pencil(513,[],1,'complex','pencil'), ...
    'unitary_hessenberg_pencil:ExplicitSize');
F=unitary_hessenberg_pencil(6,[],1,'complex','operator');
verifyError(tc,@()F(zeros(5,1),'notransp'), ...
    'um_hessenberg_pencil_apply:Dimension');
verifyError(tc,@()F(zeros(6,1),'bad'), ...
    'um_hessenberg_pencil_apply:Flag');
end
