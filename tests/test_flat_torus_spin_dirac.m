function tests=test_flat_torus_spin_dirac
%TEST_FLAT_TORUS_SPIN_DIRAC Metric-normalized 2D flat spinor baseline.
tests=functiontests(localfunctions);
end

function test_independent_sparse_dirac_exponential_and_spectrum(tc)
for nx=[3 4 5]
    for ny=[3 4]
        for beta=[0 0;pi 0;0 pi;pi pi;.32 -.17]'
            for mass=[0 .31]
                t=.63;
                W=flat_torus_spin_dirac(nx,ny,t,mass,beta,[], ...
                    'model');
                n=nx*ny;
                verifyEqual(tc,size(W.U),[2*n 2*n]);
                verifyTrue(tc,issparse(W.H));
                verifyLessThan(tc,norm(full(W.H-W.H'),'fro'), ...
                    2e-12*n);
                verifyLessThan(tc,norm(W.U-expm(-1i*t*full(W.H)), ...
                    'fro'),3e-11*n);
                verifyLessThan(tc,norm(W.U'*W.U-eye(2*n),'fro'), ...
                    3e-11*n);
                S=flat_torus_spin_dirac(nx,ny,t,mass,beta,[], ...
                    'spectral');
                verifyEqual(tc,S.finite.kind,'exact');
                lam=S.finite.eigenvalues();
                actual=eig(W.U);
                verifyLessThan(tc,max(arrayfun( ...
                    @(z)min(abs(actual-z)),lam)),2e-9);
                for k=[-3 -1 0 1 2 5]
                    verifyLessThan(tc,abs(trace(W.U^abs(k))/(2*n) ...
                         -S.finite.moments(abs(k))),3e-10*n);
                end
            end
        end
    end
end
end

function test_positive_site_metric_and_physical_M_unitarity(tc)
nx=4;ny=5;n=nx*ny;
w=.2+reshape(abs(sin(.37*(1:n))),[],1);
W=flat_torus_spin_dirac(nx,ny,.81,.26,[pi;0],w,'model');
M=W.metric;D=W.D_physical;V=W.V_physical;
verifyLessThan(tc,norm(full(D'*M-M*D),'fro'),5e-12*n);
verifyLessThan(tc,norm(full(V'*M*V-M),'fro'),5e-11*n);
verifyGreaterThan(tc,norm(V'*V-eye(2*n),'fro'),.01);
verifyLessThan(tc,norm(W.U'*W.U-eye(2*n),'fro'),2e-11*n);
U=flat_torus_spin_dirac(nx,ny,.81,.26,[pi;0],w);
V0=flat_torus_spin_dirac(nx,ny,.81,.26,[pi;0],1);
verifyLessThan(tc,norm(U-V0,'fro'),1e-13*n);
end

function test_fft_matrix_free_forward_adjoint_transpose_multirhs(tc)
for dims=[3 3;5 4;7 6;64 128]'
    nx=dims(1);ny=dims(2);n=nx*ny;
    O=flat_torus_spin_dirac(nx,ny,.73,.19,[pi;.42],[], ...
        'operator');
    verifyEqual(tc,O([],'size'),[2*n 2*n]);
    X=reshape(sin(.3*(1:6*n)),2*n,3) ...
       +1i*reshape(cos(.11*(1:6*n)),2*n,3);
    Y=O(X,'notransp');
    verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'), ...
        5e-11*norm(X,'fro'));
    verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
        5e-11*norm(X,'fro'));
    if n<=42
        U=flat_torus_spin_dirac(nx,ny,.73,.19,[pi;.42]);
        verifyLessThan(tc,norm(Y-U*X,'fro'),4e-11*norm(X,'fro'));
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            4e-11*norm(X,'fro'));
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'), ...
            4e-11*norm(X,'fro'));
        verifyGreaterThan(tc,norm(O(X,'transp') ...
            -O(X,'ctransp'),'fro'),.05);
    end
    verifyError(tc,@()O(zeros(2*n-1,1),'notransp'), ...
        'um_flat_torus_spin_dirac_apply:Dimension');
    verifyError(tc,@()O(X,'typo'), ...
        'um_flat_torus_spin_dirac_apply:Flag');
end
end

function test_spin_structures_spectral_zero_modes_and_fermion_doubling(tc)
% The four genuine torus spin structures correspond beta_x/y in {0,pi}.
for n=[5 6]
    for beta=[0 0;0 pi;pi 0;pi pi]'
        S=flat_torus_spin_dirac(n,n,1,0,beta,[], ...
            'spectral');
        z=S.finite.eigenvalues();
        verifyEqual(tc,numel(z),2*n^2);
        verifyLessThan(tc,max(abs(abs(z)-1)),2e-13);
        S2=flat_torus_spin_dirac(n,n,1,0,beta+[2*pi;0],[], ...
            'spectral');
        verifyLessThan(tc, ...
            norm(sort(angle(z))-sort(angle(S2.finite.eigenvalues()))), ...
            5e-11*n);
    end
end
% Naive centered finite differences have extra zero-energy modes on
% even-by-even grids (Nyquist sectors), not physical continuum zero modes.
S=flat_torus_spin_dirac(4,4,.63,0,[0;0],[], ...
    'spectral');
z=S.finite.eigenvalues();
verifyEqual(tc,sum(abs(z-1)<1e-12),8);
S=flat_torus_spin_dirac(5,5,.63,0,[0;0],[], ...
    'spectral');
z=S.finite.eigenvalues();
verifyEqual(tc,sum(abs(z-1)<1e-12),2);
end

function test_api_and_reproducibility(tc)
[U,p]=flat_torus_spin_dirac();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'unitary')));
S=flat_torus_spin_dirac(3,4,0,.2,[0 0],[], ...
    'spectral');
verifyEqual(tc,S.dimension,24);
verifyLessThan(tc,norm(S.finite.eigenvalues()-ones(24,1)), ...
    2e-14);
verifyError(tc,@()flat_torus_spin_dirac(2,4), ...
    'flat_torus_spin_dirac:Order');
verifyError(tc,@()flat_torus_spin_dirac(4,4,[],[],[0 0],[1 0]), ...
    'flat_torus_spin_dirac:Weights');
verifyError(tc,@()flat_torus_spin_dirac(4,4,[],[],[0]), ...
    'flat_torus_spin_dirac:Twist');
verifyError(tc,@()flat_torus_spin_dirac(20,20), ...
    'flat_torus_spin_dirac:ExplicitOrder');
verifyError(tc,@()flat_torus_spin_dirac(4,4,[],[],[],[], ...
    'wrong'),'flat_torus_spin_dirac:Mode');
end
