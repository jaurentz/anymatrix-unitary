function tests=test_cat_krawtchouk_db4
%TEST_CAT_KRAWTCHOUK_DB4 independent mathematical identities.
tests=functiontests(localfunctions);
end

function test_cat_is_unitary_and_factorizes(testCase)
for n=[1 2 3 4 5 8 9 16 19]
    for classical=[2 2;2 4;4 2]'
        a=classical(1);d=classical(2);
        U=quantum_cat_map(n,a,d);
        k=(0:n-1)';
        F=exp(-2i*pi*(k*k')/n)/sqrt(n);
        left=exp(1i*pi*d*k.^2/n);
        right=exp(1i*pi*a*k.^2/n);
        ref=bsxfun(@times,left,bsxfun(@times,F,right.'));
        verifyLessThan(testCase,norm(U-ref,'fro'),4e-11*n);
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),4e-11*n);
        verifyLessThan(testCase,max(abs(abs(U(:))-1/sqrt(n))),1e-12*n);
        Fop=quantum_cat_map(n,a,d,'operator');
        verifyEqual(testCase,Fop([],'size'),[n n]);
        X=randn(n,3)+1i*randn(n,3);
        verifyLessThan(testCase,norm(Fop(X,'notransp')-U*X,'fro'),4e-11*n);
        verifyLessThan(testCase,norm(Fop(X,'transp')-U'*X,'fro'),4e-11*n);
        verifyLessThan(testCase,norm(Fop(Fop(X,'notransp'),'transp')-X,'fro'),4e-11*n);
        S=quantum_cat_map(n,a,d,'spectral');
        verifyEqual(testCase,S.dimension,n);
        verifyEqual(testCase,S.finite.kind,'unknown');
    end
end
end

function test_cat_periodic_lift_and_pertinent_classical_matrix(testCase)
n=11;
a=2;d=4;c=a*d-1;
verifyLessThan(testCase,abs(det([a 1;c d])-1),1e-13);
verifyGreaterThan(testCase,a+d,2);
% Formula must be periodic in each integer coordinate, without a
% half-integer Bloch shift, in this admissible parity convention.
for j=-n:2*n
    for k=-n:n
        exponent=@(x,y)exp(1i*pi*(d*x.^2-2*x.*y+a*y.^2)/n);
        verifyLessThan(testCase,abs(exponent(j,k)-exponent(j+n,k)),8e-13);
        verifyLessThan(testCase,abs(exponent(j,k)-exponent(j,k+n)),8e-13);
    end
end
end

function test_cat_large_implicit_and_errors(testCase)
n=4096;
F=quantum_cat_map(n,2,2,'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-9);
verifyError(testCase,@()quantum_cat_map(8,1,3),'quantum_cat_map:Classical');
verifyError(testCase,@()quantum_cat_map(8,2,1),'quantum_cat_map:Classical');
verifyError(testCase,@()quantum_cat_map(8,2,0),'MATLAB:expectedPositive');
verifyError(testCase,@()quantum_cat_map(8,2,2,'bad'),'quantum_cat_map:Mode');
verifyError(testCase,@()F(ones(4095,1),'notransp'),'um_cat_apply:Dimension');
verifyError(testCase,@()F(ones(4096,1),'bad'),'um_cat_apply:Flag');
end

function test_krawtchouk_is_symmetric_orthogonal(testCase)
for n=[1 2 3 4 5 8 12 24 35]
    for p=[.15 .5 .83]
        U=krawtchouk_transform(n,p);
        verifyTrue(testCase,isreal(U));
        verifyLessThan(testCase,norm(U-U','fro'),2e-9*n);
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-9*n);
        verifyLessThan(testCase,norm(U^2-eye(n),'fro'),2e-9*n);
        verifyGreaterThan(testCase,min(U(1,:)),0);
        S=krawtchouk_transform(n,p,'spectral');
        verifyEqual(testCase,S.finite.kind,'exact');
        verifyEqual(testCase,numel(S.finite.eigenvalues()),n);
        verifyLessThan(testCase,abs(trace(U)-sum(S.finite.eigenvalues())),2e-8*n);
        for k=0:5
            verifyLessThan(testCase,abs(trace(U^k)/n-S.finite.moments(k)),2e-8*n);
        end
    end
end
end

function test_krawtchouk_recurrence_independent(testCase)
for p=[.27 .5 .81]
    n=17;
    U=krawtchouk_transform(n,p);
    x=(0:n-1)';
    diagv=p*(n-1-x)+(1-p)*x;
    off=-sqrt(p*(1-p)*(1:n-1)'.*(n-(1:n-1)'));
    J=diag(diagv)+diag(off,-1)+diag(off,1);
    verifyLessThan(testCase,norm(J*U-U*diag(0:n-1),'fro'),2e-10*n);
    expected=sqrt(exp(gammaln(n)-gammaln(x+1)-gammaln(n-x)) ...
        .*p.^x.*(1-p).^(n-1-x));
    verifyLessThan(testCase,norm(U(:,1)-expected),2e-10*n);
end
end

function test_krawtchouk_errors(testCase)
verifyError(testCase,@()krawtchouk_transform(7,0),'MATLAB:notGreater');
verifyError(testCase,@()krawtchouk_transform(7,1),'MATLAB:notLess');
verifyError(testCase,@()krawtchouk_transform(7,.5,'operator'),'krawtchouk_transform:Mode');
end

function test_db4_periodic_orthogonal_filterbank(testCase)
for n=[1 2 4 8 16 32 64]
    for levels=0:round(log2(n))
        U=daubechies_wavelet(n,levels);
        verifyTrue(testCase,isreal(U));
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),4e-11*n);
        F=daubechies_wavelet(n,levels,'operator');
        verifyEqual(testCase,F([],'size'),[n n]);
        X=randn(n,2)+1i*randn(n,2);
        verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),3e-11*n);
        verifyLessThan(testCase,norm(F(X,'transp')-U'*X,'fro'),3e-11*n);
        verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),3e-11*n);
        S=daubechies_wavelet(n,levels,'spectral');
        verifyEqual(testCase,S.dimension,n);
        if levels==0
            verifyEqual(testCase,S.finite.kind,'exact');
        else
            verifyEqual(testCase,S.finite.kind,'unknown');
        end
    end
end
end

function test_db4_taps_and_constant_signal(testCase)
h=[1+sqrt(3),3+sqrt(3),3-sqrt(3),1-sqrt(3)]/(4*sqrt(2));
g=[h(4),-h(3),h(2),-h(1)];
A=daubechies_wavelet(8,1);
verifyLessThan(testCase,norm(A(1,1:4)-h),1e-14);
verifyLessThan(testCase,norm(A(5,1:4)-g),1e-14);
verifyLessThan(testCase,norm(A(4,[7 8 1 2])-h),1e-14);
for n=[4 8 16 32]
    W=daubechies_wavelet(n);
    z=W*ones(n,1)/sqrt(n);
    verifyLessThan(testCase,norm(z-[1;zeros(n-1,1)]),3e-12*n);
end
end

function test_db4_large_implicit_and_errors(testCase)
n=4096;
F=daubechies_wavelet(n,12,'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),3e-9);
verifyError(testCase,@()daubechies_wavelet(6),'daubechies_wavelet:Order');
verifyError(testCase,@()daubechies_wavelet(8,4),'MATLAB:notLessEqual');
verifyError(testCase,@()daubechies_wavelet(8,[], 'bad'),'daubechies_wavelet:Mode');
verifyError(testCase,@()F(ones(n-1,1),'notransp'),'um_daubechies_apply:Dimension');
verifyError(testCase,@()F(ones(n,1),'bad'),'um_daubechies_apply:Flag');
end
