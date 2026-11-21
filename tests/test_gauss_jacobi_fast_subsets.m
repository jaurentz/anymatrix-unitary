function tests=test_gauss_jacobi_fast_subsets
%TEST_GAUSS_JACOBI_FAST_SUBSETS True FFT Gauss-Chebyshev I/II actions.
tests=functiontests(localfunctions);
end

function test_exact_dense_reference_and_fft_actions(tc)
for n=[1 2 3 4 7 8 15 32 64]
    s=(-1).^(0:n-1)';
    for a=[-.5 .5]
        Q=gauss_jacobi_transform(n,a,a);
        F=gauss_jacobi_transform(n,a,a,'operator');
        verifyTrue(tc,isa(F,'function_handle'));
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),3e-11*n);
        verifyLessThan(tc,norm(F(eye(n),'notransp')-Q,'fro'), ...
            5e-11*n);
        verifyLessThan(tc,norm(F(eye(n),'transp')-Q','fro'), ...
            5e-11*n);
        if a==-.5
            reference=orthogonal_dct(n,3)*diag(s);
        else
            reference=discrete_sine(n,1)*diag(s);
        end
        verifyLessThan(tc,norm(reference-Q,'fro'),5e-11*n);
        X=randn(n,3)+1i*randn(n,3);
        Y=randn(n,2)+1i*randn(n,2);
        verifyLessThan(tc,norm(F(X,'notransp')-Q*X,'fro'), ...
            8e-11*n);
        verifyLessThan(tc,norm(F(Y,'transp')-Q'*Y,'fro'), ...
            8e-11*n);
        verifyLessThan(tc,norm(F(Y,'ctransp')-Q'*Y,'fro'), ...
            8e-11*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'), ...
            8e-11*n);
        verifyGreaterThanOrEqual(tc,min(Q(:,1)),-3e-14);
    end
end
end

function test_large_order_actual_fft_and_adjoint(tc)
% No n-by-n matrix may be assembled by the operator path.
n=16384;
X=randn(n,4)+1i*randn(n,4);
Y=randn(n,4)+1i*randn(n,4);
for a=[-.5 .5]
    F=gauss_jacobi_transform(n,a,a,'operator');
    verifyEqual(tc,F([],'size'),[n n]);
    FX=F(X,'notransp');
    FadjY=F(Y,'transp');
    verifyLessThan(tc,norm(F(FX,'ctransp')-X,'fro'),2e-8);
    verifyLessThan(tc,abs(sum(conj(FX(:)).*Y(:))- ...
        sum(conj(X(:)).*FadjY(:))),3e-8);
end
end

function test_nonchebyshev_rejected_not_fake_fast(tc)
cases=[0 0;-.99 -.99;-.8 -.2; .25 .25; 2 2; 8 8];
for j=1:size(cases,1)
    a=cases(j,1);b=cases(j,2);
    verifyError(tc,@()gauss_jacobi_transform(12,a,b,'operator'), ...
        'gauss_jacobi_transform:FastSubset');
end
verifyError(tc,@()gauss_jacobi_transform(20,-.5,.5,'operator'), ...
    'gauss_jacobi_transform:FastSubset');
S=gauss_jacobi_transform(7,-.5,-.5,'spectral');
verifyEqual(tc,S.family,'gauss_jacobi_transform');
verifyEqual(tc,S.finite.kind,'unknown');
verifyError(tc,@()gauss_jacobi_transform(5,-.5,-.5,'bad'), ...
    'gauss_jacobi_transform:Mode');
F=gauss_jacobi_transform(6,.5,.5,'operator');
verifyError(tc,@()F(zeros(5,2),'notransp'), ...
    'um_gauss_jacobi_cheb_apply:Dimension');
verifyError(tc,@()F(zeros(6,2),'bad'), ...
    'um_gauss_jacobi_cheb_apply:Flag');
end
