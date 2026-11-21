function tests=test_periodic_rbf_propagator
%TEST_PERIODIC_RBF_PROPAGATOR Exact periodic Gaussian RBF functional calculus.
tests=functiontests(localfunctions);
end

function test_dense_expm_kernel_spectral_and_symmetry(tc)
for n=[1 2 3 4 7 11 18]
    for ell=[.25 .65 1.8]
        for t=[0 .37 -1.12]
            K=periodic_rbf_propagator(n,ell,t,'kernel');
            U=periodic_rbf_propagator(n,ell,t);
            P=periodic_rbf_propagator(n,ell,t,'parameters');
            S=periodic_rbf_propagator(n,ell,t,'spectral');
            verifyEqual(tc,size(U),[n n]);
            verifyTrue(tc,isreal(K));
            verifyLessThan(tc,norm(K-K','fro'),2e-13*n);
            verifyLessThan(tc,norm(diag(K)-ones(n,1)),2e-13*n);
            verifyLessThan(tc,norm(U-expm(1i*t*K),'fro'), ...
                3e-11*n);
            verifyLessThan(tc,norm(U'*U-eye(n),'fro'),3e-11*n);
            verifyLessThan(tc,norm(U-U.','fro'),3e-11*n);
            verifyLessThan(tc,norm(U-phase_circulant(n,t*P.kernel_eigenvalues), ...
                'fro'),3e-11*n);
            verifyLessThan(tc,norm(sort(eig(K))-sort(P.kernel_eigenvalues)), ...
                3e-11*n);
            verifyGreaterThanOrEqual(tc,min(P.kernel_eigenvalues), ...
                -3e-12*n);
            verifyLessThan(tc,abs(sum(P.kernel_eigenvalues)-n), ...
                3e-11*n);
            verifyEqual(tc,S.finite.kind,'exact');
            verifyEqual(tc,S.dimension,n);
            verifyLessThan(tc,max(arrayfun(@(z)min(abs(eig(U)-z)), ...
                S.finite.eigenvalues())),3e-9);
            for k=[-2 -1 0 1 2 5]
                if k<0,T=U'^(-k);else,T=U^k;end
                verifyLessThan(tc,abs(trace(T)/n-S.finite.moments(k)), ...
                    2e-10*n);
                verifyLessThan(tc,abs(P.moments(k)-S.finite.moments(k)), ...
                    3e-13);
            end
            if n>=3 && abs(t)>.1
                verifyGreaterThan(tc,norm(U-eye(n),'fro'),1e-3);
            end
        end
    end
end
end

function test_independent_bessel_alias_positive_fourier_spectrum(tc)
% Wrapped periodic Gaussian has strictly positive Fourier-Bessel
% coefficients: for kappa>0, K Fourier eigenvalues equal
% N sum_{q in Z} exp(-kappa) I_|m+qN|(kappa).
for n=[3 5 8 13]
    ell=.8;kappa=1/ell^2;
    P=periodic_rbf_propagator(n,ell,.43,'parameters');
    expected=zeros(n,1);
    for m=0:n-1
        q=(-12:12)';
        expected(m+1)=n*sum(besseli(abs(m+q*n),kappa,1));
    end
    verifyLessThan(tc,norm(P.kernel_eigenvalues-expected), ...
        6e-12*n);
    verifyGreaterThan(tc,min(expected),0);
    K=periodic_rbf_propagator(n,ell,.43,'kernel');
    verifyGreaterThan(tc,min(eig(K)),0);
end
end

function test_fft_operator_forward_adjoint_transpose_complex_multirhs(tc)
for n=[2 3 7 32 2048 65536]
    ell=.73;t=.43;
    O=periodic_rbf_propagator(n,ell,t,'operator');
    X=reshape(sin(.011*(1:3*n)),n,3) ...
        +1i*reshape(cos(.019*(1:3*n)),n,3);
    Y=O(X,'notransp');
    verifyEqual(tc,O([],'size'),[n n]);
    verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'), ...
        3e-11*norm(X,'fro'));
    verifyLessThan(tc,norm(O(Y,'transp')-O(Y,'notransp'),'fro'), ...
        3e-11*norm(X,'fro'));
    verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
        3e-11*norm(X,'fro'));
    if n<=32
        U=periodic_rbf_propagator(n,ell,t);
        verifyLessThan(tc,norm(Y-U*X,'fro'),3e-11*norm(X,'fro'));
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            3e-11*norm(X,'fro'));
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'), ...
            3e-11*norm(X,'fro'));
        if n>2,verifyGreaterThan(tc, ...
                norm(O(X,'transp')-O(X,'ctransp'),'fro'),.001);end
    end
    verifyError(tc,@()O(zeros(n+1,1),'notransp'), ...
        'um_periodic_rbf_apply:Dimension');
    verifyError(tc,@()O(X,'bad'), ...
        'um_periodic_rbf_apply:Flag');
end
end

function test_validation_and_special_limits(tc)
[U,props]=periodic_rbf_propagator();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyError(tc,@()periodic_rbf_propagator(0), ...
    'MATLAB:expectedPositive');
verifyError(tc,@()periodic_rbf_propagator(4,-.2), ...
    'MATLAB:expectedPositive');
verifyError(tc,@()periodic_rbf_propagator(4,.4,[], ...
    'bad'),'periodic_rbf_propagator:Mode');
verifyError(tc,@()periodic_rbf_propagator(2049,.4,.2), ...
    'periodic_rbf_propagator:DenseOrder');
O=periodic_rbf_propagator(2049,.4,.2,'operator');
verifyEqual(tc,O([],'size'),[2049 2049]);
P=periodic_rbf_propagator(8,.5,0,'parameters');
verifyLessThan(tc,norm(P.propagator_eigenvalues-ones(8,1)), ...
    3e-13);
S=periodic_rbf_propagator(8,.5,0,'spectral');
verifyEqual(tc,S.finite.moments([-4 -1 0 1 5]),ones(1,5));
verifyError(tc,@()S.finite.moments(1.5),'MATLAB:expectedInteger');
end
