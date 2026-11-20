function tests=test_cat_krawtchouk_daubechies
%TEST_CAT_KRAWTCHOUK_DAUBECHIES Independent source-based identities.
tests=functiontests(localfunctions);
end

function test_cat_hannay_berry_formula_all_parities(testCase)
for n=[1 2 3 4 5 7 8 12 17]
    for s=[1 2 4]
        Q=quantum_cat(n,s);
        verifyEqual(testCase,size(Q),[n n]);
        verifyLessThan(testCase,norm(Q'*Q-eye(n),'fro'),3e-11*n);
        verifyLessThan(testCase,max(abs(abs(Q(:))-1/sqrt(n))),3e-12);
        j=(0:n-1)';
        E=exp(1i*pi/4)*exp(2i*pi*(j.^2*ones(1,n) ...
          -j*j' +s*ones(n,1)*(j'.^2))/n)/sqrt(n);
        verifyLessThan(testCase,norm(Q-E,'fro'),3e-11*n);
        F=quantum_cat(n,s,'operator');
        verifyEqual(testCase,F([],'size'),[n n]);
        X=randn(n,3)+1i*randn(n,3);
        verifyLessThan(testCase,norm(F(X,'notransp')-Q*X,'fro'),2e-11*n);
        verifyLessThan(testCase,norm(F(X,'transp')-Q'*X,'fro'),2e-11*n);
        verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),2e-11*n);
        S=quantum_cat(n,s,'spectral');
        verifyEqual(testCase,S.finite.kind,'unknown');
        verifyEqual(testCase,S.dimension,n);
    end
end
end

function test_cat_classical_symplectic_constraints(testCase)
for s=[1 2 4 8]
    T=[2*s 1;4*s-1 2];
    verifyEqual(testCase,round(det(T)),1);
    verifyGreaterThan(testCase,trace(T),2);
    verifyEqual(testCase,mod(T(1,1),2),0);
    verifyEqual(testCase,mod(T(2,2),2),0);
    verifyEqual(testCase,mod(T(1,2),2),1);
    verifyEqual(testCase,mod(T(2,1),2),1);
end
end

function test_cat_fft_large(testCase)
n=8192;
F=quantum_cat(n,3,'operator');
X=randn(n,3)+1i*randn(n,3);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-9);
end

function test_cat_validation(testCase)
verifyError(testCase,@()quantum_cat(5,0),'quantum_cat:BadS');
end

function test_krawtchouk_tridiagonal_diagonalization(testCase)
for n=[1 2 3 6 12 24]
    for p=[.2 .5 .8]
        Q=krawtchouk_transform(n,p);
        verifyTrue(testCase,isreal(Q));
        verifyEqual(testCase,size(Q),[n n]);
        verifyLessThan(testCase,norm(Q'*Q-eye(n),'fro'),5e-11*n);
        N=n-1;
        if n==1
            verifyEqual(testCase,Q,1);
            continue
        end
        j=(0:N)';
        diagv=N*p+j*(1-2*p);
        off=sqrt(p*(1-p)*(1:N)'.*(N-(1:N)'+1));
        J=diag(diagv)+diag(off,1)+diag(off,-1);
        x=0:N;
        verifyLessThan(testCase,norm(J*Q-Q*diag(x),'fro'),2e-10*n);
        verifyGreaterThanOrEqual(testCase,min(Q(1,:)),-1e-13);
        w=arrayfun(@(k)nchoosek(N,k)*p^k*(1-p)^(N-k),x);
        verifyLessThan(testCase,norm(Q(1,:).^2-w),2e-11*n);
        S=krawtchouk_transform(n,p,'spectral');
        verifyEqual(testCase,S.finite.kind,'unknown');
        verifyEqual(testCase,S.dimension,n);
    end
end
end

function test_krawtchouk_distinction_from_jacobi_spectrum(testCase)
Q=krawtchouk_transform(9,.4);
z=eig(Q);
verifyLessThan(testCase,max(abs(abs(z)-1)),2e-10);
verifyGreaterThan(testCase,norm(sort(real(z))-(0:8)'),1);
end

function test_krawtchouk_validation(testCase)
verifyError(testCase,@()krawtchouk_transform(8,0),'MATLAB:expectedScalar');
verifyError(testCase,@()krawtchouk_transform(8,1),'MATLAB:expectedScalar');
verifyError(testCase,@()krawtchouk_transform(8,.5,'operator'), ...
    'krawtchouk_transform:Mode');
end

function test_daubechies_orthogonality_and_operator(testCase)
for order=2:4
    for n=[1 2 4 8 16 32 64]
        for levels=0:round(log2(n))
            W=daubechies_wavelet(n,order,levels);
            verifyTrue(testCase,isreal(W));
            verifyEqual(testCase,size(W),[n n]);
            verifyLessThan(testCase,norm(W'*W-eye(n),'fro'),5e-11*n);
            F=daubechies_wavelet(n,order,levels,'operator');
            verifyEqual(testCase,F([],'size'),[n n]);
            X=randn(n,3)+1i*randn(n,3);
            verifyLessThan(testCase,norm(F(X,'notransp')-W*X,'fro'),2e-11*n);
            verifyLessThan(testCase,norm(F(X,'transp')-W'*X,'fro'),2e-11*n);
            verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),2e-11*n);
        end
    end
end
end

function test_daubechies_periodic_filter_and_dc_annihilation(testCase)
n=32;
for order=2:4
    W=daubechies_wavelet(n,order);
    verifyLessThan(testCase,norm(W*ones(n,1)-[sqrt(n);zeros(n-1,1)]),3e-11*n);
    verifyGreaterThan(testCase,norm(W-haar_wavelet(n),'fro'),.1);
    S=daubechies_wavelet(n,order,[],'spectral');
    verifyEqual(testCase,S.finite.kind,'unknown');
end
S0=daubechies_wavelet(n,2,0,'spectral');
verifyEqual(testCase,S0.finite.kind,'exact');
verifyEqual(testCase,numel(S0.finite.eigenvalues()),n);
verifyEqual(testCase,S0.finite.moments(-2:3),ones(1,6));
end

function test_daubechies_implicit_large(testCase)
n=8192;
F=daubechies_wavelet(n,4,[],'operator');
X=randn(n,3)+1i*randn(n,3);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),1e-9);
end

function test_daubechies_validation(testCase)
verifyError(testCase,@()daubechies_wavelet(7),'daubechies_wavelet:Order');
verifyError(testCase,@()daubechies_wavelet(8,1),'MATLAB:expectedScalar');
verifyError(testCase,@()daubechies_wavelet(8,5),'MATLAB:expectedScalar');
verifyError(testCase,@()daubechies_wavelet(8,2,4),'MATLAB:expectedScalar');
verifyError(testCase,@()daubechies_wavelet(8,2,2,'invalid'),'daubechies_wavelet:Mode');
F=daubechies_wavelet(8,2,[],'operator');
verifyError(testCase,@()F(ones(7,1),'notransp'),'um_daubechies_apply:Dimension');
verifyError(testCase,@()F(ones(8,1),'invalid'),'um_daubechies_apply:Flag');
end
