function tests=test_gauss_jacobi_fast_subsets
%TEST_GAUSS_JACOBI_FAST_SUBSETS True FFT Gauss-Chebyshev I--IV actions.
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

function test_third_fourth_kind_direct_formula_and_dense_reference(tc)
% Independent trigonometric Gauss formulas, not a DCT/DST alias.
% The n quadrature nodes are in ASCENDING x order (descending theta).
for n=[1 2 3 4 7 8 15 32 64 96]
    M=2*n+1;
    j=(n:-1:1)';
    k=0:n-1;
    for pair=[-.5 .5; .5 -.5]'
        a=pair(1);b=pair(2);
        if a<0
            theta=(2*j-1)*pi/M;
            reference=(2/sqrt(M))*cos(theta*(k+.5));
        else
            theta=2*j*pi/M;
            reference=(2/sqrt(M))*sin(theta*(k+.5));
        end
        Q=gauss_jacobi_transform(n,a,b);
        F=gauss_jacobi_transform(n,a,b,'operator');
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(Q-reference,'fro'),8e-11*n);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),3e-11*n);
        verifyLessThan(tc,norm(F(eye(n),'notransp')-Q,'fro'), ...
            8e-11*n);
        verifyLessThan(tc,norm(F(eye(n),'transp')-Q','fro'), ...
            8e-11*n);
        X=complex(reshape(sin((1:3*n)*.34),n,3), ...
                  reshape(cos((1:3*n)*.27),n,3));
        Y=complex(reshape(cos((1:2*n)*.41),n,2), ...
                  reshape(sin((1:2*n)*.29),n,2));
        verifyLessThan(tc,norm(F(X,'notransp')-Q*X,'fro'), ...
            2e-10*n);
        verifyLessThan(tc,norm(F(Y,'transp')-Q'*Y,'fro'), ...
            2e-10*n);
        verifyLessThan(tc,norm(F(Y,'ctransp')-Q'*Y,'fro'), ...
            2e-10*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X, ...
            'fro'),2e-10*n);
        verifyLessThan(tc,norm(imag(F(real(X),'notransp')),'fro'), ...
            2e-12*n);
        verifyGreaterThanOrEqual(tc,min(Q(:,1)),-3e-13);
    end
end
end

function test_large_order_chebyshev_third_fourth_fft(tc)
% The 16384-by-16384 evaluation matrix is NEVER assembled.
n=16384;
v=(1:n)';
X=[sin(.013*v)+1i*cos(.017*v),cos(.021*v)+1i*sin(.033*v)];
Y=[cos(.025*v)+1i*sin(.019*v),sin(.031*v)+1i*cos(.023*v)];
for pair=[-.5 .5;.5 -.5]'
    F=gauss_jacobi_transform(n,pair(1),pair(2),'operator');
    FX=F(X,'notransp');
    FY=F(Y,'ctransp');
    verifyEqual(tc,size(FX),size(X));
    verifyLessThan(tc,norm(F(FX,'ctransp')-X,'fro'),3e-8);
    verifyLessThan(tc,abs(sum(conj(FX(:)).*Y(:)) ...
        -sum(conj(X(:)).*FY(:))),3e-8);
end
end

function test_nonchebyshev_rejected_not_fake_fast(tc)
cases=[0 0;-.99 -.99;-.8 -.2; .25 .25; 2 2; 8 8];
for j=1:size(cases,1)
    a=cases(j,1);b=cases(j,2);
    verifyError(tc,@()gauss_jacobi_transform(12,a,b,'operator'), ...
        'gauss_jacobi_transform:FastSubset');
end
F3=gauss_jacobi_transform(20,-.5,.5,'operator');
F4=gauss_jacobi_transform(20,.5,-.5,'operator');
verifyEqual(tc,F3([],'size'),[20 20]);
verifyEqual(tc,F4([],'size'),[20 20]);
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
