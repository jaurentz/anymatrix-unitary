function tests=test_nonuniform_fem_wave
%TEST_NONUNIFORM_FEM_WAVE Sparse noncommuting FE wave / Gauss geometry.
tests=functiontests(localfunctions);
end

function test_independent_element_assembly_and_noncommutator(tc)
for n=[2 3 4 8 12]
    for omega0=[0 .8]
        coeff=nonuniform_fem_wave(n,2,1.3,omega0,.02, ...
            'midpoint','coefficients');
        [Mass,K,widths,speeds]=independent_assembly(n,2,1.3,omega0);
        verifyLessThan(tc,norm(full(coeff.mass-Mass),'fro'),2e-13);
        verifyLessThan(tc,norm(full(coeff.stiffness-K),'fro'),2e-11);
        verifyLessThan(tc,norm(coeff.nodes-[0;cumsum(widths)]),1e-13);
        verifyLessThan(tc,norm(coeff.element_speeds-speeds),1e-14);
        verifyGreaterThan(tc,min(eig(full(coeff.mass))),0);
        verifyGreaterThan(tc,min(eig(full(coeff.stiffness))),0);
        verifyGreaterThan(tc,norm(full(coeff.mass*coeff.stiffness ...
            -coeff.stiffness*coeff.mass),'fro'),1e-5);
        verifyGreaterThan(tc,coeff.commutator_fro,1e-5);
        Q=nonuniform_fem_wave(n,2,1.3,omega0,.02, ...
            'midpoint','quadratic');
        verifyLessThan(tc,norm(full(Q.A0-Q.stiffness),'fro'),1e-13);
        verifyLessThan(tc,norm(full(Q.A2-Q.mass),'fro'),1e-13);
        verifyEqual(tc,nnz(Q.A1),0);
    end
end
end

function test_Gauss_physical_metric_symplecticity_and_orthogonality(tc)
for n=[2 3 6 9]
    for method={'midpoint','gauss2','gauss3'}
        name=method{1};h=.07/(n+1);
        P=nonuniform_fem_wave(n,1,1,.7,h,name,'pencil');
        U=nonuniform_fem_wave(n,1,1,.7,h,name);
        J=[zeros(n),eye(n);-eye(n),zeros(n)];
        verifyTrue(tc,isreal(U));
        verifyEqual(tc,size(U),[2*n 2*n]);
        verifyLessThan(tc,norm(U'*U-eye(2*n),'fro'), ...
            3e-9*n);
        verifyLessThan(tc,norm(P.U-U,'fro'),3e-12*n);
        verifyLessThan(tc,norm(P.Mstep'*P.G*P.Mstep-P.G,'fro')/ ...
            norm(P.G,'fro'),3e-9);
        verifyLessThan(tc,norm(P.Mstep'*J*P.Mstep-J,'fro'), ...
            3e-9*n);
        verifyLessThan(tc,norm(P.T'*P.T-P.G,'fro')/ ...
            norm(P.G,'fro'),3e-11);
        verifyGreaterThan(tc,min(eig(P.G)),0);
        verifyLessThan(tc,norm(P.T*P.Mstep/P.T-U,'fro'), ...
            3e-10*n);
        verifyLessThan(tc,norm(P.A'* (P.B\P.A)-P.B,'fro')/ ...
            norm(P.B,'fro'),3e-9);
        verifyLessThan(tc,norm(U'*P.normalizedSymplectic*U ...
             -P.normalizedSymplectic,'fro'),3e-9*n);
        verifyLessThan(tc,norm(P.physicalGenerator'*P.G ...
              +P.G*P.physicalGenerator,'fro'),1e-10*n);
        verifyLessThan(tc,norm(P.physicalGenerator'*J ...
              +J*P.physicalGenerator,'fro'),1e-10*n);
        verifyGreaterThan(tc,norm(U'*J*U-J,'fro'),.02);
        lambda=eig(P.A,P.B);
        verifyLessThan(tc,max(abs(abs(lambda)-1)),3e-8);
    end
end
end

function test_sparse_rational_action_and_real_orthogonal_adjoint(tc)
for n=[2 4 7 15]
    for method={'midpoint','gauss2','gauss3'}
        name=method{1};h=.06/(n+1);
        O=nonuniform_fem_wave(n,[],[],[],h,name,'operator');
        U=nonuniform_fem_wave(n,[],[],[],h,name);
        X=reshape(sin(.21*(1:6*n)),2*n,3) ...
           +1i*reshape(cos(.33*(1:6*n)),2*n,3);
        verifyEqual(tc,O([],'size'),[2*n 2*n]);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'), ...
            1e-8*n);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            1e-8*n);
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'), ...
            1e-8*n);
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
            1e-8*n);
        R=O(real(X),'notransp');
        verifyTrue(tc,isreal(R));
        verifyLessThan(tc,norm(R-U*real(X),'fro'),1e-8*n);
    end
end
end

function test_large_nonuniform_sparse_LU_actions(tc)
for n=[160 384]
    O=nonuniform_fem_wave(n,1.1,1.3,.5,[], ...
        'gauss3','operator');
    x=complex(cos((1:2*n)'*.31),sin((1:2*n)'*.43));
    y=O(x,'notransp');
    verifyEqual(tc,O([],'size'),[2*n 2*n]);
    verifyLessThan(tc,norm(O(y,'ctransp')-x), ...
        3e-8*norm(x));
    verifyLessThan(tc,abs(norm(y)-norm(x)), ...
        3e-8*norm(x));
    verifyError(tc,@()O(ones(2*n-1,1),'notransp'), ...
        'um_nonuniform_fem_apply:Dimension');
    verifyError(tc,@()O(x,'nope'), ...
        'um_nonuniform_fem_apply:Flag');
end
end

function test_independent_Pade_and_expm_physical_step(tc)
n=5;L=1.3;c=1.2;omega0=.6;
coeff=nonuniform_fem_wave(n,L,c,omega0,.01, ...
    'midpoint','coefficients');
M=full(coeff.mass);K=full(coeff.stiffness);
A0=[zeros(n),M\eye(n);-K,zeros(n)];
h=.09/(n+1);
for method={'midpoint','gauss2','gauss3'}
    name=method{1};
    P=nonuniform_fem_wave(n,L,c,omega0,h,name,'pencil');
    coeffs=independent_pade(name);
    Z=h*A0;
    F=zeros(2*n);B=zeros(2*n);
    for j=1:numel(coeffs)
        F=F+coeffs(j)*Z^(j-1);
        B=B+coeffs(j)*(-Z)^(j-1);
    end
    expected=B\F;
    verifyLessThan(tc,norm(P.Mstep-expected,'fro'),2e-9);
    verifyLessThan(tc,norm(P.physicalGenerator-A0,'fro'),1e-11);
    e1=norm(P.Mstep-expm(h*A0),'fro');
    Phalf=nonuniform_fem_wave(n,L,c,omega0,h/2,name,'pencil');
    e2=norm(Phalf.Mstep-expm(.5*h*A0),'fro');
    verifyGreaterThan(tc,e1,4*e2);
end
end

function test_homogeneous_uniform_limit_exact_finite_spectra(tc)
uniform=struct('mesh_amplitude',0,'speed_amplitude',0);
for n=[2 4 8]
    for method={'midpoint','gauss2','gauss3'}
        name=method{1};
        h=.035/(n+1);
        S=nonuniform_fem_wave(n,1,1,.8,h,name, ...
            'spectral',uniform);
        R=geometric_fem_wave(n,1,1,.8,h,name,'spectral');
        C=nonuniform_fem_wave(n,1,1,.8,h,name, ...
            'coefficients',uniform);
        verifyLessThan(tc,norm(full(C.mass* ...
            C.stiffness-C.stiffness*C.mass),'fro'),2e-12);
        verifyEqual(tc,S.finite.kind,'exact');
        verifyLessThan(tc,norm(S.finite.eigenvalues() ...
           -R.finite.eigenvalues()),4e-12*n);
        verifyLessThan(tc,norm(S.frequencies ...
           -R.frequencies),2e-12*n);
        U=nonuniform_fem_wave(n,1,1,.8,h,name, ...
            'matrix',uniform);
        specEig=S.finite.eigenvalues();
        Ueig=eig(U);
        for j=1:numel(specEig)
            verifyLessThan(tc,min(abs(Ueig-specEig(j))),2e-8);
        end
        for k=-3:4
            verifyLessThan(tc,abs(trace(U^abs(k))/(2*n) ...
                -S.finite.moments(abs(k))),3e-10);
        end
    end
end
unknown=nonuniform_fem_wave(5,1,1,.8,.03, ...
    'midpoint','spectral');
verifyEqual(tc,unknown.finite.kind,'unknown');
verifyEmpty(tc,unknown.finite.eigenvalues);
verifyEmpty(tc,unknown.finite.moments);
verifyGreaterThan(tc,unknown.mass_commutator_fro,1e-5);
end

function test_explicit_heterogeneous_mesh_option_and_conditioning(tc)
n=5;L=1.1;
w=[.09;.2;.15;.11;.3;.25];w=w*L/sum(w);
ce=[1.2;1.1;.7;1.4;1.3;.9];
o=struct('element_lengths',w,'element_speeds',ce);
C=nonuniform_fem_wave(n,L,1,.7,.012, ...
    'midpoint','coefficients',o);
verifyLessThan(tc,norm(C.element_lengths-w),1e-13);
verifyLessThan(tc,norm(C.element_speeds-ce),1e-13);
verifyGreaterThan(tc,C.commutator_fro,.001);
O=nonuniform_fem_wave(n,L,1,.7,.012, ...
    'midpoint','operator',o);
U=nonuniform_fem_wave(n,L,1,.7,.012, ...
    'midpoint','matrix',o);
verifyLessThan(tc,norm(O(eye(2*n),'notransp')-U,'fro'),2e-10);
w=[1e-16;1e-16;.55;1e-16;.45-3e-16];
ill=struct('element_lengths',w,'element_speeds',ones(5,1));
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02, ...
    'midpoint','coefficients',ill), ...
    'nonuniform_fem_wave:Conditioning');
end

function test_properties_and_errors(tc)
[Z,p]=nonuniform_fem_wave();
verifyEmpty(tc,Z);
verifyTrue(tc,any(strcmp(p,'orthogonal')));
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyEqual(tc,size(nonuniform_fem_wave(4)),[8 8]);
verifyError(tc,@()nonuniform_fem_wave(1),'nonuniform_fem_wave:Order');
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02,'bad'), ...
    'nonuniform_fem_wave:Method');
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02, ...
    'midpoint','bad'),'nonuniform_fem_wave:Mode');
verifyError(tc,@()nonuniform_fem_wave(129), ...
    'nonuniform_fem_wave:ExplicitOrder');
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02, ...
    'midpoint','matrix',struct('mesh_amplitude',1)), ...
    'nonuniform_fem_wave:Options');
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02, ...
    'midpoint','matrix',struct('unknown',1)), ...
    'nonuniform_fem_wave:Options');
verifyError(tc,@()nonuniform_fem_wave(4,1,1,.5,.02, ...
    'midpoint','matrix',struct('element_lengths',ones(5,1))), ...
    'nonuniform_fem_wave:Mesh');
end

function [M,K,widths,speeds]=independent_assembly(n,L,c,omega0)
e=(1:n+1)';
w=1+.27*sin(sqrt(2)*e+.4);
widths=L*w/sum(w);
speeds=c*(1+.34*cos(.81*e+.2));
M=zeros(n);K=zeros(n);
for j=1:n+1
    de=widths(j);v=speeds(j);
    localM=de/6*[2,1;1,2];
    localK=v*v/de*[1,-1;-1,1];
    nodes=[j-1,j];
    for a=1:2
        for b=1:2
            if nodes(a)>=1 && nodes(a)<=n ...
                    &&nodes(b)>=1 &&nodes(b)<=n
                M(nodes(a),nodes(b))= ...
                    M(nodes(a),nodes(b))+localM(a,b);
                K(nodes(a),nodes(b))= ...
                    K(nodes(a),nodes(b))+localK(a,b);
            end
        end
    end
end
K=K+omega0^2*M;
M=sparse(M);K=sparse(K);
end

function coeff=independent_pade(method)
switch method
    case 'midpoint',coeff=[1 .5];
    case 'gauss2',coeff=[1 .5 1/12];
    case 'gauss3',coeff=[1 .5 .1 1/120];
end
end
