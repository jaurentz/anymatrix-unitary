function tests=test_daubechies_orders
%TEST_DAUBECHIES_ORDERS Extension of the canonical periodic wavelet family.
tests=functiontests(localfunctions);
end

function test_orders_orthogonality_and_implicit(tc)
for order=2:4
    for n=[1 2 4 8 16 32 64]
        for levels=0:round(log2(n))
            W=daubechies_wavelet(n,levels,'matrix',order);
            verifyTrue(tc,isreal(W));
            verifyLessThan(tc,norm(W'*W-eye(n),'fro'),5e-11*n);
            F=daubechies_wavelet(n,levels,'operator',order);
            X=randn(n,3)+1i*randn(n,3);
            verifyLessThan(tc,norm(F(X,'notransp')-W*X,'fro'),2e-11*n);
            verifyLessThan(tc,norm(F(X,'transp')-W'*X,'fro'),2e-11*n);
            verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),2e-11*n);
            if order==2
                verifyLessThan(tc,norm(W-daubechies_wavelet(n,levels),'fro'),1e-13*n);
            end
        end
    end
end
end

function test_dc_and_legacy_signature(tc)
n=64;
h=[1+sqrt(3),3+sqrt(3),3-sqrt(3),1-sqrt(3)]/(4*sqrt(2));
W=daubechies_wavelet(8,1);
verifyLessThan(tc,norm(W(1,1:4)-h),1e-13);
for order=2:4
    U=daubechies_wavelet(n,[],[],order);
    verifyLessThan(tc,norm(U*ones(n,1)-[sqrt(n);zeros(n-1,1)]),4e-11*n);
    S=daubechies_wavelet(n,[],'spectral',order);
    verifyEqual(tc,S.finite.kind,'unknown');
end
S0=daubechies_wavelet(n,0,'spectral',4);
verifyEqual(tc,S0.finite.kind,'exact');
end

function test_large_db4_and_errors(tc)
n=8192;
F=daubechies_wavelet(n,[],'operator',4);
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'),1e-9);
verifyError(tc,@()daubechies_wavelet(8,1,'matrix',1),'daubechies_wavelet:Filter');
verifyError(tc,@()daubechies_wavelet(8,1,'matrix',5),'daubechies_wavelet:Filter');
verifyError(tc,@()daubechies_wavelet(8,4,'matrix',3),'daubechies_wavelet:Levels');
end
