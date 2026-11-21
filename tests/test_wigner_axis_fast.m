function tests=test_wigner_axis_fast
%TEST_WIGNER_AXIS_FAST Exact sparse special-angle spin rotations (#41).
tests=functiontests(localfunctions);
end

function test_source_convention_against_independent_matrix_exponential(tc)
for ell=[0 .5 1 1.5 2 3.5 5]
    n=round(2*ell)+1;
    m=(-ell:ell)';
    if ell>0
        c=sqrt((ell-m(1:end-1)).*(ell+m(1:end-1)+1));
        K=(diag(c,1)-diag(c,-1))/2;
    else,K=0;end
    for q=-4:5
        beta=q*pi;alpha=.37;gamma=-.54;
        reference=diag(exp(-1i*alpha*m))*expm(beta*K) ...
            *diag(exp(-1i*gamma*m));
        U=wigner_rotation(ell,alpha,beta,gamma);
        verifyLessThan(tc,norm(full(U)-reference,'fro'), ...
            2e-12*n);
        verifyLessThan(tc,norm(full(U)'*full(U)-eye(n),'fro'), ...
            2e-12*n);
        O=wigner_rotation(ell,alpha,beta,gamma,'operator');
        X=reshape(sin(.13*(1:3*n)),n,3) ...
            +1i*reshape(cos(.27*(1:3*n)),n,3);
        verifyEqual(tc,O([],'size'),[n n]);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'), ...
            4e-12*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            4e-12*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'), ...
            4e-12*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X, ...
            'fro'),4e-12*max(1,norm(X,'fro')));
        S=wigner_rotation(ell,alpha,beta,gamma,'spectral');
        lam=S.finite.eigenvalues();
        actual=eig(full(U));
        verifyLessThan(tc,max(arrayfun(@(z)min(abs(actual-z)), ...
            lam)),3e-10);
    end
end
end

function test_SU2_double_cover_and_exact_diagonal_antidiagonal(tc)
for ell=[1 2 3.5 4.5]
    n=round(2*ell)+1;
    eta=(-1)^round(2*ell);
    U0=wigner_rotation(ell,0,0,0);
    U2=wigner_rotation(ell,0,2*pi,0);
    U4=wigner_rotation(ell,0,4*pi,0);
    verifyLessThan(tc,norm(full(U0-eye(n)),'fro'),2e-13);
    verifyLessThan(tc,norm(full(U2-eta*eye(n)),'fro'),2e-13);
    verifyLessThan(tc,norm(full(U4-eye(n)),'fro'),2e-13);
    P=wigner_rotation(ell,0,pi,0);
    verifyLessThan(tc,norm(full(P*P-eta*eye(n)),'fro'),2e-13);
end
for ell=[.5 1 3 4.5]
    n=round(2*ell)+1;
    alpha=.6;gamma=.19;
    D=wigner_rotation(ell,alpha,0,gamma);
    verifyLessThan(tc,norm(D-diag(diag(D)),'fro'),2e-13);
    P=wigner_rotation(ell,alpha,pi,gamma);
    mask=fliplr(eye(n));
    verifyLessThan(tc,norm(full(P).* (1-mask),'fro'),2e-13);
end
end

function test_large_spin_linear_storage_and_actions(tc)
for ell=[512 2048.5]
    n=round(2*ell)+1;
    for beta=[0 pi 2*pi 3*pi]
        O=wigner_rotation(ell,.27,beta,-.36,'operator');
        X=reshape(sin(.013*(1:3*n)),n,3) ...
            +1i*reshape(cos(.019*(1:3*n)),n,3);
        Y=O(X,'notransp');
        verifyEqual(tc,O([],'size'),[n n]);
        verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'), ...
            2e-12*norm(X,'fro'));
        verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
            2e-12*norm(X,'fro'));
        U=wigner_rotation(ell,.27,beta,-.36);
        verifyTrue(tc,issparse(U));
        verifyEqual(tc,nnz(U),n);
        verifyLessThan(tc,norm(U*X-Y,'fro'),2e-12*norm(X,'fro'));
    end
end
end

function test_non_axis_rejected_without_approximation(tc)
verifyError(tc,@()wigner_rotation(2,.17,.3,-.2,'operator'), ...
    'wigner_rotation:Operator');
verifyError(tc,@()wigner_rotation(2,.17,pi+1e-8,-.2,'operator'), ...
    'wigner_rotation:Operator');
verifyError(tc,@()wigner_rotation(3,.17,[],[],'wrong'), ...
    'wigner_rotation:Mode');
O=wigner_rotation(3,.2,pi,.4,'operator');
verifyError(tc,@()O(zeros(6,1),'notransp'), ...
    'um_wigner_axis_apply:Dimension');
verifyError(tc,@()O(ones(7,1),'wrong'), ...
    'um_wigner_axis_apply:Flag');
end
