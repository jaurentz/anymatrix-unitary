function tests=test_spherical_cap_commuting
%TEST_SPHERICAL_CAP_COMMUTING Source-matched tridiagonal Slepian blocks.
tests=functiontests(localfunctions);
end

function test_commutation_spectral_order_and_orthogonality(tc)
for L=[1 2 4 7 12 20]
    n=L^2;
    for cap=[.21 .93 2.32]
        C=spherical_slepian(L,cap,'commuting');
        P=spherical_slepian(L,cap,'parameters');
        K=P.concentration;
        O=spherical_slepian(L,cap,'operator');
        Q=O(eye(n),'ctransp');
        lam=C.concentration_eigenvalues;
        verifyEqual(tc,O([],'size'),[n n]);
        verifyEqual(tc,numel(C.tridiagonal_blocks),L);
        verifyEqual(tc,numel(C.blocks),2*L-1);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'), ...
            2e-10*n);
        verifyLessThan(tc,norm(Q'*K*Q-diag(lam),'fro'), ...
            3e-9*n);
        verifyLessThan(tc,norm(lam ...
            -P.concentration_eigenvalues),6e-9*n);
        verifyLessThan(tc,abs(sum(lam) ...
            -n*(1-cos(cap))/2),1e-9*n);
        verifyGreaterThanOrEqual(tc,min(lam),-1e-10);
        verifyLessThanOrEqual(tc,max(lam),1+1e-10);
        for m=0:L-1
            ell=(m:L-1)';
            h=numel(ell);
            T=full(C.tridiagonal_blocks{m+1});
            verifyEqual(tc,size(T),[h h]);
            verifyLessThan(tc,norm(T-T','fro'),1e-11);
            verifyLessThan(tc,norm(T-triu(tril(T,1),-1),'fro'), ...
                1e-13);
            % Published Grünbaum matrix differs by an overall sign.
            verifyLessThan(tc,norm(diag(T) ...
                -cos(cap)*ell.*(ell+1)),1e-11);
            if h>1
                k=(m:L-2)';
                off=(L^2-(k+1).^2).*sqrt( ...
                    ((k+1).^2-m^2)./((2*k+1).*(2*k+3)));
                verifyLessThan(tc,norm(diag(T,1)-off),1e-11);
            end
            ids=find(P.azimuthal_order==m & P.harmonic_kind==(m>0));
            if m==0,ids=find(P.azimuthal_order==0);end
            verifyEqual(tc,numel(ids),h);
            Kblock=K(ids,ids);
            verifyLessThan(tc,norm(T*Kblock-Kblock*T,'fro'), ...
                1e-9*max(1,norm(T,'fro')));
            % Interlaced cosine/sine copies are identical operators.
            if m>0
                other=find(P.azimuthal_order==m & P.harmonic_kind==2);
                verifyLessThan(tc,norm(K(other,other)-Kblock,'fro'),1e-11);
            end
        end
    end
end
end

function test_forward_adjoint_transpose_and_multirhs(tc)
for L=[3 9 16 48]
    n=L^2;cap=.72;
    O=spherical_slepian(L,cap,'operator');
    X=reshape(sin(.018*(1:3*n)),n,3) ...
        +1i*reshape(cos(.023*(1:3*n)),n,3);
    Z=O(X,'notransp');
    verifyLessThan(tc,norm(O(Z,'ctransp')-X,'fro'), ...
        8e-11*norm(X,'fro'));
    verifyLessThan(tc,norm(O(Z,'transp')-X,'fro'), ...
        8e-11*norm(X,'fro'));
    verifyLessThan(tc,abs(norm(Z,'fro')-norm(X,'fro')), ...
        8e-11*norm(X,'fro'));
    if L<=16
        U=spherical_slepian(L,cap,'matrix');
        K=spherical_slepian(L,cap,'concentration');
        lam=spherical_slepian(L,cap,'commuting');
        verifyLessThan(tc,norm((O(eye(n),'ctransp'))' ...
             *K*O(eye(n),'ctransp') ...
             -diag(lam.concentration_eigenvalues),'fro'), ...
             2e-9*n);
        % Nondegenerate/mildly clustered eigenvectors have matching
        % squared entry magnitudes, not necessarily their global signs.
        A=O(eye(n),'notransp');
        for j=1:n
            % Compare eigenvector projectors only at distinguishable
            % concentration factors; avoid arbitrary near-unit doublets.
            if min(abs(lam.concentration_eigenvalues(j) ...
                    -lam.concentration_eigenvalues(setdiff(1:n,j))))>1e-8
                verifyLessThan(tc,norm(abs(A(j,:)) ...
                    -abs(U(j,:))),3e-7);
            end
        end
    end
    verifyError(tc,@()O(zeros(n-1,1),'notransp'), ...
        'um_spherical_cap_slepian_apply:Dimension');
    verifyError(tc,@()O(X,'wrong'), ...
        'um_spherical_cap_slepian_apply:Flag');
end
end

function test_full_empty_caps_and_reproducibility(tc)
for L=[1 2 8 24]
    n=L^2;
    for cap=[0 pi]
        O=spherical_slepian(L,cap,'operator');
        C=spherical_slepian(L,cap,'commuting');
        X=reshape(cos((1:2*n)*.014),n,2);
        verifyLessThan(tc,norm(O(X,'notransp')-X,'fro'), ...
            3e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(X,'ctransp')-X,'fro'), ...
            3e-11*max(1,norm(X,'fro')));
        verifyEqual(tc,C.concentration_eigenvalues, ...
            double(cap==pi)*ones(n,1));
    end
end
A=spherical_slepian(13,.83,'commuting');
B=spherical_slepian(13,.83,'commuting');
verifyEqual(tc,A.concentration_eigenvalues, ...
    B.concentration_eigenvalues,'AbsTol',0);
for m=1:13
    verifyEqual(tc,A.tridiagonal_blocks{m}, ...
        B.tridiagonal_blocks{m},'AbsTol',0);
end
verifyError(tc,@()spherical_slepian(25,.6,'matrix'), ...
    'spherical_slepian:DenseOrder');
verifyError(tc,@()spherical_slepian(65,.6,'operator'), ...
    'MATLAB:notLessEqual');
S=spherical_slepian(48,.71,'spectral');
verifyEqual(tc,S.dimension,48^2);
verifyEqual(tc,S.finite.kind,'unknown');
end
