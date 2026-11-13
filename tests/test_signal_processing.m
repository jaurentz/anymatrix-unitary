function tests=test_signal_processing
%TEST_SIGNAL_PROCESSING Independent identities for finite DSP unitaries.
tests=functiontests(localfunctions);
end

function test_dct_and_hartley(testCase)
for n=[1 2 3 4 7 16]
    for kind=1:4
        if n==1 && kind==1,continue;end
        Q=orthogonal_dct(n,kind);
        verifyTrue(testCase,isreal(Q));
        unitary(testCase,Q,3e-12*n);
        if any(kind==[1 4])
            verifyLessThan(testCase,norm(Q-Q','fro'),3e-12*n);
            verifyLessThan(testCase,norm(Q*Q-eye(n),'fro'),3e-12*n);
        end
    end
    H=discrete_hartley(n);
    unitary(testCase,H,3e-12*n);
    verifyLessThan(testCase,norm(H-H','fro'),3e-12*n);
    verifyLessThan(testCase,norm(H*H-eye(n),'fro'),3e-12*n);
end
verifyLessThan(testCase,norm(orthogonal_dct(7,3)-orthogonal_dct(7,2)','fro'),1e-12);
verifyError(testCase,@()orthogonal_dct(1,1),'orthogonal_dct:Order');
end

function test_cazac_autocorrelation_and_phases(testCase)
for n=[1 2 3 4 5 7 8 13]
    for r=1:min(n,5)
        if gcd(r,n)~=1,continue;end
        U=cazac_circulant(n,r);
        unitary(testCase,U,2e-11*n);
        verifyLessThan(testCase,max(abs(abs(U(:))-1/sqrt(n))),3e-12);
        lambda=fft(U(:,1));
        for k=1:min(n,4)
            verifyLessThan(testCase,abs(trace(U^k)-sum(lambda.^k)),2e-10*n);
        end
    end
end
verifyError(testCase,@()cazac_circulant(6,2),'cazac_circulant:Coprime');
end

function test_fresnel_fourier_reference(testCase)
for n=[1 2 3 4 5 8 11 16]
    U=fresnel_transform(n);
    unitary(testCase,U,2e-11*n);
    c=U(:,1);lambda=fft(c);
    for k=1:min(n,4)
        verifyLessThan(testCase,abs(trace(U^k)-sum(lambda.^k)),2e-10*n);
    end
    for j=1:n
        verifyLessThan(testCase,norm(U(:,j)-circshift(c,j-1)),5e-12*n);
    end
end
end

function test_affine_and_zak(testCase)
for n=[1 2 5 8 12]
    k=(0:n-1)';
    F=exp(-2i*pi*(k*k')/n)/sqrt(n);
    verifyLessThan(testCase,norm(affine_fourier(n,0,0)-F,'fro'),2e-12*n);
    unitary(testCase,affine_fourier(n,.031,-.07),2e-11*n);
    for p=1:n
        if mod(n,p)~=0,continue;end
        Z=finite_zak(n,p);
        unitary(testCase,Z,2e-11*n);
        q=n/p;
        t=(0:q-1)';
        Fq=exp(-2i*pi*(t*t')/q)/sqrt(q);
        verifyLessThan(testCase,norm(Z-kron(Fq,eye(p)),'fro'),1e-12*n);
    end
end
verifyError(testCase,@()finite_zak(7,2),'finite_zak:Factor');
end

function test_weyl_shifts(testCase)
for n=[1 2 3 7 10]
    for p=[-2 0 1 3]
        for q=[-3 0 1 4]
            U=time_frequency_shift(n,p,q);
            verifyTrue(testCase,issparse(U));
            verifyEqual(testCase,nnz(U),n);
            unitary(testCase,U,2e-12*n);
            j=0:n-1; row=mod(j+q,n);
            expected=exp(2*pi*1i*mod(p,n)*row/n);
            actual=full(U(sub2ind([n,n],row+1,j+1)));
            verifyLessThan(testCase,norm(actual(:)-expected(:)),2e-12*n);
        end
    end
end
end

function test_periodic_hilbert_completion(testCase)
for n=[1 2 3 4 5 8 15]
    for dc=[-1 1]
        for nyq=[-1 1]
            U=periodic_hilbert(n,dc,nyq);
            verifyTrue(testCase,isreal(U));
            unitary(testCase,U,3e-11*n);
            d=zeros(n,1); d(1)=dc;
            d(2:ceil(n/2))=-1i;
            d(floor(n/2)+2:end)=1i;
            if mod(n,2)==0,d(n/2+1)=nyq;end
            for k=1:min(n,5)
                verifyLessThan(testCase,abs(trace(U^k)-sum(d.^k)),3e-10*n);
            end
        end
    end
end
verifyError(testCase,@()periodic_hilbert(4,0),'periodic_hilbert:Signs');
end

function unitary(testCase,U,tol)
n=size(U,1);
verifyEqual(testCase,size(U),[n n]);
verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),tol);
end
