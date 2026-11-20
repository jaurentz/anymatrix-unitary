function tests=test_paraunitary_fir
tests=functiontests(localfunctions);
end

function test_frequency_polyphase_and_winding(tc)
for d=[1 2 3 4]
    m=3;
    V=randn(d,m)+1i*randn(d,m);
    for j=1:m,V(:,j)=V(:,j)/norm(V(:,j));end
    U0=discrete_fourier(d);
    C=paraunitary_fir(7,d,V,U0,[],'polyphase');
    verifyEqual(tc,size(C.taps),[d d m+1]);
    verifyEqual(tc,C.winding,-m);
    for omega=linspace(-pi,pi,19)
        E=paraunitary_fir(7,d,V,U0,omega,'frequency');
        verifyLessThan(tc,norm(E'*E-eye(d),'fro'),4e-11*d);
        Z=zeros(d);
        for j=0:m,Z=Z+C.taps(:,:,j+1)*exp(-1i*j*omega);end
        verifyLessThan(tc,norm(Z-E,'fro'),4e-11*d);
        verifyLessThan(tc,abs(det(E)*exp(1i*m*omega)-det(U0)), ...
            4e-11*d);
    end
    for lag=-m:m
        M=zeros(d);
        for r=0:m
            s=r+lag;
            if s>=0&&s<=m
                M=M+C.taps(:,:,r+1)'*C.taps(:,:,s+1);
            end
        end
        expected=zeros(d);
        if lag==0,expected=eye(d);end
        verifyLessThan(tc,norm(M-expected,'fro'),4e-11*d);
    end
end
end

function test_periodic_operator_and_small_fourier_spectra(tc)
for L=[1 2 3 4 8]
    for d=[1 2 3]
        V=randn(d,2)+1i*randn(d,2);
        U0=discrete_fourier(d);
        U=paraunitary_fir(L,d,V,U0);
        verifyEqual(tc,size(U),[L*d L*d]);
        verifyLessThan(tc,norm(U'*U-eye(L*d),'fro'),3e-11*L*d);
        F=paraunitary_fir(L,d,V,U0,[],'operator');
        X=randn(L*d,3)+1i*randn(L*d,3);
        verifyEqual(tc,F([],'size'),[L*d L*d]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),3e-11*L*d);
        verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),3e-11*L*d);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),3e-11*L*d);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            3e-11*L*d);
        % Independently Fourier-transform the periodic channel block.
        T=kron(discrete_fourier(L),eye(d));
        W=T*U*T';
        E=zeros(d*L);
        for k=0:L-1
            pos=k*d+(1:d);
            E(pos,pos)=paraunitary_fir(L,d,V,U0,2*pi*k/L,'frequency');
        end
        verifyLessThan(tc,norm(W-E,'fro'),3e-11*L*d);
        S=paraunitary_fir(L,d,V,U0,[],'spectral');
        verifyEqual(tc,S.finite.kind,'exact');
        verifyEqual(tc,S.dimension,d*L);
        z=S.finite.eigenvalues();
        verifyEqual(tc,numel(z),d*L);
        for e=eig(U).'
            verifyLessThan(tc,min(abs(z-e)),3e-9);
        end
        for k=-3:4
            verifyLessThan(tc,abs(S.finite.moments(k)- ...
                trace(U^k)/(L*d)),5e-10);
        end
    end
end
end

function test_large_matrix_free_and_errors(tc)
L=2048;d=3;
F=paraunitary_fir(L,d,[],[],[],'operator');
X=randn(L*d,2)+1i*randn(L*d,2);
verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'),4e-10);
verifyError(tc,@()paraunitary_fir(2,2,zeros(2,2)), ...
    'paraunitary_fir:Directions');
verifyError(tc,@()paraunitary_fir(2,2,ones(3,2)), ...
    'paraunitary_fir:Directions');
verifyError(tc,@()paraunitary_fir(2,2,[],2*eye(2)), ...
    'paraunitary_fir:U0');
verifyError(tc,@()paraunitary_fir(2,2,[],[],Inf), ...
    'MATLAB:expectedFinite');
verifyError(tc,@()paraunitary_fir(2,2,[],[],[],'invalid'), ...
    'paraunitary_fir:Mode');
verifyError(tc,@()F(ones(L*d-1,1),'notransp'), ...
    'um_paraunitary_apply:Dimension');
verifyError(tc,@()F(ones(L*d,1),'bad'), ...
    'um_paraunitary_apply:Flag');
end
