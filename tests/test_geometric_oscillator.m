function tests=test_geometric_oscillator
%TEST_GEOMETRIC_OSCILLATOR Spectral and quadratic-invariant verification.
tests=functiontests(localfunctions);
end

function test_orthogonal_symplectic_and_exact_phase_pairs(tc)
for n=[1,2,5,10]
    for method={'midpoint','gauss2','gauss3','verlet'}
        name=method{1};
        U=geometric_oscillator(n,0.45,1.1,0.7,name);
        ref=geometric_oscillator(n,0.45,1.1,0.7,name,'spectral');
        verifyEqual(tc,size(U),[2*n,2*n]);
        verifyEqual(tc,ref.dimension,2*n);
        verifyEqual(tc,ref.kind,'exact');
        verifyLessThan(tc,norm(U'*U-eye(2*n),'fro'),2e-11*n);
        J=[zeros(n),eye(n);-eye(n),zeros(n)];
        verifyLessThan(tc,norm(U'*J*U-J,'fro'),2e-11*n);
        eigs=eig(U);
        for j=1:numel(ref.eigenvalues)
            verifyLessThan(tc,min(abs(eigs-ref.eigenvalues(j))),2e-8);
        end
    end
end
end

function test_positive_metric_pencil_and_physical_step(tc)
n=7;h=0.5;om=1.3;kappa=0.8;
for method={'midpoint','gauss2','gauss3','verlet'}
    m=method{1};
    P=geometric_oscillator(n,h,om,kappa,m,'pencil');
    U=geometric_oscillator(n,h,om,kappa,m);
    G=P.G;M=P.M;T=P.T;
    verifyEqual(tc,P.classification,'circle_spectrum_only');
    verifyLessThan(tc,norm(P.U-U,'fro'),1e-12);
    verifyGreaterThan(tc,min(eig(G)),0);
    verifyLessThan(tc,norm(M'*G*M-G,'fro'),2e-9*norm(G,'fro'));
    verifyLessThan(tc,norm(P.A-G*M,'fro'),1e-12*norm(G,'fro'));
    verifyLessThan(tc,norm(P.B-G,'fro'),1e-12);
    verifyLessThan(tc,norm(P.A'/P.B*P.A-P.B,'fro'), ...
        2e-9*norm(G,'fro'));
    verifyLessThan(tc,norm(T*M/T-U,'fro'),2e-9);
    J=[zeros(n),eye(n);-eye(n),zeros(n)];
    verifyLessThan(tc,norm(M'*J*M-J,'fro'),2e-9);
    roots=eig(P.A,P.B);
    ref=eig(U);
    for j=1:numel(ref)
        verifyLessThan(tc,min(abs(roots-ref(j))),2e-8);
    end
end
end

function test_verlet_palindromic_recurrence(tc)
for n=[1 2 6]
    h=0.4;om=1.0;kappa=0.6;
    p=geometric_oscillator(n,h,om,kappa,'verlet','polynomial');
    verifyEqual(tc,p.A0,eye(n),'AbsTol',0);
    verifyEqual(tc,p.A2,p.A0,'AbsTol',0);
    verifyLessThan(tc,norm(p.A1-p.A1','fro'),1e-12);
    C=[zeros(n),eye(n);-p.A0,-p.A1];
    eigs=eig(C);
    z=geometric_oscillator(n,h,om,kappa,'verlet','spectral');
    for j=1:2*n
        verifyLessThan(tc,min(abs(eigs-z.eigenvalues(j))),2e-8);
    end
end
end

function test_gauss_scalar_pade_and_verlet_frequency(tc)
n=8;h=0.4;om=1.2;kappa=0.9;
for s=1:3
    methods={'midpoint','gauss2','gauss3'};
    ref=geometric_oscillator(n,h,om,kappa,methods{s},'spectral');
    w=ref.frequencies;
    switch s
        case 1, a=[1 1/2];
        case 2, a=[1 1/2 1/12];
        case 3, a=[1 1/2 1/10 1/120];
    end
    z=polyval(fliplr(a),-1i*h*w)./polyval(fliplr(a),1i*h*w);
    verifyLessThan(tc,max(abs(ref.eigenvalues(1:n)-z)),1e-13);
end
ref=geometric_oscillator(n,h,om,kappa,'verlet','spectral');
verifyLessThan(tc,max(abs(ref.eigenvalues(1:n)- ...
    exp(-2i*asin(h*ref.frequencies/2)))),1e-13);
end

function test_defaults_and_validation(tc)
[u,props]=geometric_oscillator();
verifyEmpty(tc,u);
verifyTrue(tc,any(strcmp(props,'orthogonal')));
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyEqual(tc,size(geometric_oscillator(4)),[8 8]);
verifyError(tc,@()geometric_oscillator(4,1.5,1,1,'verlet'), ...
    'geometric_oscillator:Stability');
verifyError(tc,@()geometric_oscillator(4,.5,1,1,'gauss2','polynomial'), ...
    'geometric_oscillator:Polynomial');
verifyError(tc,@()geometric_oscillator(4,.5,1,1,'bad'), ...
    'geometric_oscillator:Method');
end

function test_matrix_free_fft_actions(tc)
for n=[1,2,7,13]
    for method={'midpoint','gauss2','gauss3','verlet'}
        name=method{1};
        U=geometric_oscillator(n,.42,1.2,.5,name);
        O=geometric_oscillator(n,.42,1.2,.5,name,'operator');
        verifyTrue(tc,isa(O,'function_handle'));
        verifyEqual(tc,O([],'size'),[2*n,2*n]);
        X=randn(2*n,3)+1i*randn(2*n,3);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'), ...
            3e-11*max(1,n));
        verifyLessThan(tc,norm(O(X,'transp')-U'*X,'fro'), ...
            3e-11*max(1,n));
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            3e-11*max(1,n));
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
            3e-11*max(1,n));
        verifyLessThan(tc,norm(O(speye(2*n),'notransp')-U,'fro'), ...
            3e-11*max(1,n));
    end
end
end

function test_large_matrix_free(tc)
n=2048;
O=geometric_oscillator(n,.42,1.2,.5,'gauss3','operator');
X=randn(2*n,2)+1i*randn(2*n,2);
Y=O(X,'notransp');
verifyEqual(tc,O([],'size'),[2*n,2*n]);
verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'),2e-10);
verifyError(tc,@()O(ones(2*n-1,1),'notransp'), ...
    'um_geometric_oscillator_apply:Dimension');
verifyError(tc,@()O(X,'invalid'), ...
    'um_geometric_oscillator_apply:Flag');
end
