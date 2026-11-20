function tests=test_wavelet_packets
tests=functiontests(localfunctions);
end

function test_full_trees_periodic_aliases_and_adjoint(tc)
for order=2:4
    for n=[1 2 4 8 16 32]
        for levels=0:round(log2(n))
            W=daubechies_wavelet(n,levels,'matrix',order,'full');
            verifyEqual(tc,size(W),[n n]);
            verifyTrue(tc,isreal(W));
            verifyLessThan(tc,norm(W'*W-eye(n),'fro'),5e-11*n);
            F=daubechies_wavelet(n,levels,'operator',order,'full');
            X=randn(n,3)+1i*randn(n,3);
            verifyEqual(tc,F([],'size'),[n n]);
            verifyLessThan(tc,norm(F(X,'notransp')-W*X,'fro'),5e-11*n);
            verifyLessThan(tc,norm(F(X,'transp')-W'*X,'fro'),5e-11*n);
            verifyLessThan(tc,norm(F(X,'ctransp')-W'*X,'fro'),5e-11*n);
            verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),5e-11*n);
        end
    end
end
end

function test_legacy_and_incomplete_trees(tc)
for order=2:4
    for n=[8 16 32]
        depth=round(log2(n));
        A=daubechies_wavelet(n,depth,'matrix',order);
        B=daubechies_wavelet(n,depth,'matrix',order,2.^(0:depth-1));
        verifyLessThan(tc,norm(A-B,'fro'),2e-11*n);
        W=daubechies_wavelet(n,3,'matrix',order,[1 2 3 6]);
        verifyLessThan(tc,norm(W'*W-eye(n),'fro'),5e-11*n);
        F=daubechies_wavelet(n,3,'operator',order,[1 2 3 6]);
        X=randn(n,3)+1i*randn(n,3);
        verifyLessThan(tc,norm(F(X,'notransp')-W*X,'fro'),5e-11*n);
        verifyLessThan(tc,norm(F(X,'transp')-W'*X,'fro'),5e-11*n);
    end
end
S=daubechies_wavelet(16,3,'spectral',4,'none');
verifyEqual(tc,S.finite.kind,'exact');
verifyLessThan(tc,norm(daubechies_wavelet(16,3,'matrix',4,'none')-eye(16),'fro'),1e-14);
end

function test_custom_orthonormal_quadrature_mirror_filter(tc)
h=[1;1]/sqrt(2);
for n=[2 4 8 32]
    depth=round(log2(n));
    W=daubechies_wavelet(n,depth,'matrix',h,'full');
    verifyLessThan(tc,norm(W'*W-eye(n),'fro'),4e-12*n);
    F=daubechies_wavelet(n,depth,'operator',h,'full');
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(tc,norm(F(F(X,'notransp'),'transp')-X,'fro'),4e-12*n);
end
U=daubechies_wavelet(16,3,'matrix',2,'full');
V=daubechies_wavelet(16,3,'matrix',4,'full');
verifyGreaterThan(tc,norm(U-V,'fro'),.1);
S=daubechies_wavelet(16,3,'spectral',4,'full');
verifyEqual(tc,S.finite.kind,'unknown');
end

function test_packet_and_filter_errors(tc)
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',2,[1 4]), ...
    'daubechies_wavelet:Packet');
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',2,[1 2 2]), ...
    'daubechies_wavelet:Packet');
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',2,[1 2 8]), ...
    'daubechies_wavelet:Packet');
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',2,'bogus'), ...
    'daubechies_wavelet:Packet');
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',[1 0],'full'), ...
    'daubechies_wavelet:Filter');
verifyError(tc,@()daubechies_wavelet(16,3,'matrix',ones(4,1)/2,'full'), ...
    'daubechies_wavelet:Filter');
F=daubechies_wavelet(16,3,'operator',2,'full');
verifyError(tc,@()F(ones(15,2),'notransp'), ...
    'um_wavelet_packet_apply:Dimension');
verifyError(tc,@()F(ones(16,2),'invalid'), ...
    'um_wavelet_packet_apply:Flag');
end
