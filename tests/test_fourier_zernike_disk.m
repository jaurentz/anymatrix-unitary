function tests=test_fourier_zernike_disk
%TEST_FOURIER_ZERNIKE_DISK Exact finite radial Jacobi/angular Fourier map.
tests=functiontests(localfunctions);
end

function test_unitary_matrix_and_forward_adjoint_transpose(tc)
for nr=[1 2 3 5 7]
    for np=[2 3 4 5 8]
        N=nr*np;
        U=fourier_zernike_disk(nr,np);
        F=fourier_zernike_disk(nr,np,'operator');
        verifyEqual(tc,size(U),[N N]);
        verifyEqual(tc,F([],'size'),[N N]);
        verifyLessThan(tc,norm(U'*U-eye(N),'fro'), ...
            7e-12*max(1,N));
        X=reshape(sin(.23*(1:N*3)),N,3) ...
            +1i*reshape(cos(.17*(1:N*3)),N,3);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'), ...
            1e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'), ...
            1e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(X,'transp')-U.'*X,'fro'), ...
            1e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X, ...
            'fro'),1e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(F(X,'ctransp'),'notransp')-X, ...
            'fro'),1e-11*max(1,norm(X,'fro')));
    end
end
end

function test_independent_physical_zernike_jacobi_orthonormal_samples(tc)
% Independent finite binomial sum for classical
% R_(m+2k)^m(r)=r^m P_k^(0,m)(2*r^2-1).
for nr=[2 3 5 8]
    for np=[3 4 7]
        M=fourier_zernike_disk(nr,np,'model');
        U=fourier_zernike_disk(nr,np);
        for col=1:np
            m=abs(M.signed_modes(col));
            r=M.radial_nodes_by_mode(:,col);
            w=M.physical_radial_weights_by_mode(:,col);
            verifyGreaterThan(tc,min(w),0);
            verifyGreaterThan(tc,min(r),0);
            verifyLessThan(tc,max(r),1);
            Q=zeros(nr);
            for k=0:nr-1
                p=independent_jacobi0m(k,m,2*r.^2-1);
                % L2(r dr) radial Zernike functions:
                % normalization sqrt(2*(2k+m+1)).
                physical=sqrt(2*(2*k+m+1))*r.^m.*p;
                Q(:,k+1)=sqrt(w).*physical;
            end
            verifyLessThan(tc,norm(M.radial_blocks{col}-Q,'fro'), ...
                8e-11*nr);
            verifyLessThan(tc,norm(Q'*Q-eye(nr),'fro'), ...
                8e-11*nr);
            verifyEqual(tc,M.radial_degree_by_mode(:,col), ...
                (m+2*(0:nr-1))');
            for k=1:nr
                C=zeros(nr,np);C(k,col)=1;
                X=C(:);
                Y=reshape(U*X,nr,np);
                T=Q(:,k)*exp(1i*M.signed_modes(col) ...
                    *M.angular_angles)/sqrt(np);
                verifyLessThan(tc,norm(Y-T,'fro'), ...
                    1e-10*max(1,norm(T,'fro')));
            end
        end
        verifyFalse(tc,M.is_common_physical_grid);
        if np>=3
            verifyGreaterThan(tc,norm(M.radial_nodes_by_mode(:,1) ...
                -M.radial_nodes_by_mode(:,2)),1e-4);
        end
    end
end
end

function test_large_implicit_transform_no_dense_unitary(tc)
nr=48;np=64;N=nr*np;
F=fourier_zernike_disk(nr,np,'operator');
X=reshape(sin(.021*(1:2*N)),N,2) ...
    +1i*reshape(cos(.013*(1:2*N)),N,2);
Y=F(X,'notransp');
verifyEqual(tc,F([],'size'),[N N]);
verifyLessThan(tc,norm(F(Y,'ctransp')-X,'fro'), ...
    4e-9*max(1,norm(X,'fro')));
verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
    4e-9*max(1,norm(X,'fro')));
verifyError(tc,@()fourier_zernike_disk(48,64,'matrix'), ...
    'fourier_zernike_disk:DenseOrder');
end

function test_spectral_unknown_and_explicit_validation(tc)
[U,props]=fourier_zernike_disk();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(props,'unitary')));
S=fourier_zernike_disk(4,5,'spectral');
verifyEqual(tc,S.dimension,20);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
verifyError(tc,@()fourier_zernike_disk(0,5), ...
    'fourier_zernike_disk:RadialOrder');
verifyError(tc,@()fourier_zernike_disk(4,1), ...
    'fourier_zernike_disk:AngularOrder');
verifyError(tc,@()fourier_zernike_disk(90,64), ...
    'fourier_zernike_disk:Order');
verifyError(tc,@()fourier_zernike_disk(4,5,'wrong'), ...
    'fourier_zernike_disk:Mode');
F=fourier_zernike_disk(4,5,'operator');
verifyError(tc,@()F(zeros(19,1),'notransp'), ...
    'um_fourier_zernike_apply:Dimension');
verifyError(tc,@()F(zeros(20,1),'wrong'), ...
    'um_fourier_zernike_apply:Flag');
end

function p=independent_jacobi0m(k,m,x)
% Polynomial from finite binomial formula, not Golub-Welsch eigenvectors.
p=zeros(size(x));
for s=0:k
    p=p+nchoosek(k,s)*nchoosek(k+m,s) ...
        *((x-1)/2).^s.*((x+1)/2).^(k-s);
end
end
