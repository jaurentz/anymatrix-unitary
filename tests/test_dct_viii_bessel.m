function tests=test_dct_viii_bessel
%TEST_DCT_VIII_BESSEL Finite exact nu=-1/2 Hankel / odd-period DCT-VIII.
tests=functiontests(localfunctions);
end

function test_bessel_kernel_is_exact_dct_viii(testCase)
for n=[1 2 3 5 8 13 32]
    Q=orthogonal_dct(n,8);
    a=pi*((1:n)'-0.5); % positive zeros of J_{-1/2}
    c=(n+0.5)*pi;
    d=abs(besselj(0.5,a));
    H=(2/c)*besselj(-0.5,(a*a')/c)./(d*d');
    verifyLessThan(testCase,norm(H-Q,'fro'),2e-11*n);
    verifyLessThan(testCase,norm(Q'*Q-eye(n),'fro'),3e-12*n);
    verifyLessThan(testCase,norm(Q-Q','fro'),3e-12*n);
    verifyLessThan(testCase,norm(Q*Q-eye(n),'fro'),3e-12*n);
    verifyTrue(testCase,isreal(Q));
    if n>1
        verifyGreaterThan(testCase,norm(Q-orthogonal_dct(n,4),'fro'),1e-3);
    end
    S=orthogonal_dct(n,8,'spectral');
    z=S.finite.eigenvalues();
    verifyLessThan(testCase,abs(trace(Q)-sum(z)),5e-12*n);
end
end

function test_fft_actions_and_independent_cosine(testCase)
for n=[1 2 3 5 6 15 32]
    Q=orthogonal_dct(n,8);
    k=(0:n-1)';
    R=(2/sqrt(2*n+1))*cos(2*pi*(k+0.5).*(k'+0.5)/(2*n+1));
    verifyLessThan(testCase,norm(Q-R,'fro'),1e-12*n);
    F=orthogonal_dct(n,8,'operator');
    verifyEqual(testCase,F([],'size'),[n n]);
    X=randn(n,3)+1i*randn(n,3);
    for flag={'notransp','transp','ctransp'}
        verifyLessThan(testCase,norm(F(X,flag{1})-Q*X,'fro'),3e-11*n);
    end
    verifyLessThan(testCase,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),4e-11*n);
    verifyLessThan(testCase,norm(F(eye(n),'notransp')-Q,'fro'),4e-11*n);
end
end

function test_large_fft_and_validation(testCase)
n=8192;
F=orthogonal_dct(n,8,'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-9);
verifyError(testCase,@()orthogonal_dct(5,5),'MATLAB:notLessEqual');
verifyError(testCase,@()orthogonal_dct(5,6),'MATLAB:notLessEqual');
verifyError(testCase,@()orthogonal_dct(5,7),'MATLAB:notLessEqual');
verifyError(testCase,@()F(ones(n-1,1),'notransp'),'um_dct_apply:Dimension');
verifyError(testCase,@()F(ones(n,1),'invalid'),'um_dct_apply:Flag');
S=orthogonal_dct(8,8,'spectral');
verifyEqual(testCase,S.dimension,8);
verifyEqual(testCase,S.finite.kind,'exact');
z=S.finite.eigenvalues();
verifyEqual(testCase,sum(abs(z-1)<eps),4);
verifyEqual(testCase,sum(abs(z+1)<eps),4);
verifyEqual(testCase,S.finite.atoms.weights,[.5;.5],'AbsTol',eps);
verifyEqual(testCase,S.finite.moments(0:5),[1 0 1 0 1 0],'AbsTol',eps);
Sodd=orthogonal_dct(7,8,'spectral');
verifyEqual(testCase,Sodd.finite.atoms.weights,[4/7;3/7],'AbsTol',eps);
verifyEqual(testCase,Sodd.finite.moments(0:3),[1 1/7 1 1/7],'AbsTol',eps);
verifyEqual(testCase,orthogonal_dct(8,4,'spectral').finite.kind,'unknown');
end
