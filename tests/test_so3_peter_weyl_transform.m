function tests=test_so3_peter_weyl_transform
%TEST_SO3_PETER_WEYL_TRANSFORM Haar normalization and Wigner FFT sampling.
tests=functiontests(localfunctions);
end

function test_fixed_wigner_frequency_square_block(tc)
for L=[1 2 3 5 8 12]
    for m=-min(L-1,3):min(L-1,3)
        for n=-min(L-1,3):min(L-1,3)
            P=L-max(abs(m),abs(n));
            Q=so3_peter_weyl_transform(L,m,n);
            verifyEqual(tc,size(Q),[P P]);
            verifyTrue(tc,isreal(Q));
            verifyLessThan(tc,norm(Q'*Q-eye(P),'fro'),2e-10*P);
            a=abs(m-n);b=abs(m+n);
            ref=(-1)^max(0,m-n)*gauss_jacobi_transform(P,a,b);
            verifyLessThan(tc,norm(Q-ref,'fro'),2e-12*P);
            S=so3_peter_weyl_transform(L,m,n,'spectral');
            verifyEqual(tc,S.finite.kind,'unknown');
            verifyEqual(tc,S.dimension,P);
            par=so3_peter_weyl_transform(L,m,n,'parameters');
            verifyEqual(tc,par.fixed_block_dimension,P);
            verifyEqual(tc,par.sampling_rows,L*(2*L-1)^2);
            verifyEqual(tc,par.sampling_columns,L*(4*L^2-1)/3);
            verifyLessThan(tc,abs(sum(par.beta_weights_dx)-2),3e-12*L);
            verifyLessThan(tc,abs(sum(par.normalized_haar_weights) ...
                *(2*L-1)^2-1),3e-12*L);
        end
    end
end
end

function test_peter_weyl_full_haar_isometry(tc)
for L=1:5
    A=so3_peter_weyl_transform(L,0,0,'sampling');
    F=so3_peter_weyl_transform(L,0,0,'sampling_operator');
    C=L*(4*L^2-1)/3;R=L*(2*L-1)^2;
    verifyEqual(tc,size(A),[R C]);
    verifyEqual(tc,F([],'size'),[R C]);
    verifyLessThan(tc,norm(A'*A-eye(C),'fro'), ...
        1e-10*max(C,1));
    X=randn(C,3)+1i*randn(C,3);
    Y=randn(R,2)+1i*randn(R,2);
    verifyLessThan(tc,norm(F(X,'notransp')-A*X,'fro'), ...
        1e-10*max(C,1));
    verifyLessThan(tc,norm(F(Y,'transp')-A'*Y,'fro'), ...
        1e-10*max(C,1));
    verifyLessThan(tc,norm(F(Y,'ctransp')-A'*Y,'fro'), ...
        1e-10*max(C,1));
    verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
        1e-10*max(C,1));
    verifyLessThan(tc,abs(norm(F(X,'notransp'),'fro')-norm(X,'fro')), ...
        1e-10*max(C,1));
    if L>1,verifyGreaterThan(tc,R,C);end
end
end

function test_wigner_matrix_element_and_haar_constant(tc)
for L=[1 2 3 4]
    P=so3_peter_weyl_transform(L,0,0,'parameters');
    A=so3_peter_weyl_transform(L,0,0,'sampling');
    M=2*L-1;
    expected=repmat(sqrt(P.beta_weights_dx/(2*M^2)),M^2,1);
    verifyLessThan(tc,norm(A(:,1)-expected),4e-12*L);
    if L<3,continue;end
    ell=2;m=1;n=-1;
    col=0;
    for k=0:ell-1,col=col+(2*k+1)^2;end
    col=col+(m+ell)*(2*ell+1)+(n+ell)+1;
    for ib=[1 L]
        for ia=[1 M]
            ig=M;
            idx=ib+(ia-1)*L+(ig-1)*L*M;
            aa=P.alpha_nodes(ia);
            bb=acos(P.beta_cos_nodes(ib));
            gg=P.gamma_nodes(ig);
            D=wigner_rotation(ell,aa,bb,gg);
            val=sqrt(P.beta_weights_dx(ib)/(2*M^2))* ...
                sqrt(2*ell+1)*D(ell+1+m,ell+1+n);
            verifyLessThan(tc,abs(A(idx,col)-val),4e-11*L);
        end
    end
end
end

function test_so3_z_axis_left_right_translation(tc)
for L=[2 3 5]
    M=2*L-1;
    F=so3_peter_weyl_transform(L,0,0,'sampling_operator');
    C=L*(4*L^2-1)/3;
    x=randn(C,2)+1i*randn(C,2);
    left=x;right=x;ids=0;
    for ell=0:L-1
        for m=-ell:ell
            for n=-ell:ell
                ids=ids+1;
                left(ids,:)=exp(-2i*pi*m/M)*x(ids,:);
                right(ids,:)=exp(-2i*pi*n/M)*x(ids,:);
            end
        end
    end
    y=F(x,'notransp');
    yl=F(left,'notransp');
    yr=F(right,'notransp');
    for rhs=1:2
        V=reshape(y(:,rhs),L,M,M);
        expectedL=circshift(V,[0 -1 0]);
        expectedR=circshift(V,[0 0 -1]);
        verifyLessThan(tc,norm(yl(:,rhs)-expectedL(:)), ...
            7e-11*C);
        verifyLessThan(tc,norm(yr(:,rhs)-expectedR(:)), ...
            7e-11*C);
    end
end
end

function test_large_operator_and_invalid_contracts(tc)
L=12;F=so3_peter_weyl_transform(L,1,-2,'sampling_operator');
sz=F([],'size');
verifyEqual(tc,sz,[L*(2*L-1)^2,L*(4*L^2-1)/3]);
X=randn(sz(2),2)+1i*randn(sz(2),2);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
    2e-9*L);
verifyError(tc,@()so3_peter_weyl_transform(8,0,0,'sampling'), ...
    'so3_peter_weyl_transform:SamplingSize');
verifyError(tc,@()so3_peter_weyl_transform(3,3,0), ...
    'so3_peter_weyl_transform:Order');
verifyError(tc,@()so3_peter_weyl_transform(3,0,0,'bad'), ...
    'so3_peter_weyl_transform:Mode');
verifyError(tc,@()F(zeros(sz(2)-1,2),'notransp'), ...
    'um_so3_sampling_apply:Dimension');
verifyError(tc,@()F(zeros(sz(1),1),'bad'), ...
    'um_so3_sampling_apply:Flag');
end
