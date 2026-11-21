function tests=test_finite_koopman_torus
%TEST_FINITE_KOOPMAN_TORUS Exact modular Koopman unitary permutation.
tests=functiontests(localfunctions);
end

function test_independent_affine_map_and_sparse_orthogonality(tc)
for m=[2 3 4 5 8]
    n=m*m;
    cases={struct('M',[2 1;1 1],'b',[0;0]), ...
      struct('M',[1 1;0 1],'b',[0;0]), ...
      struct('M',eye(2),'b',[1;0]), ...
      struct('M',[0 -1;1 0],'b',[1;1])};
    for t=1:numel(cases)
        M=cases{t}.M;shift=cases{t}.b;
        U=finite_koopman_torus(m,M,'matrix',shift);
        d=finite_koopman_torus(m,M,'map',shift);
        verifyEqual(tc,size(U),[n,n]);
        verifyTrue(tc,issparse(U));
        verifyTrue(tc,isreal(U));
        verifyLessThan(tc,norm(U'*U-speye(n),'fro'),2e-14);
        verifyEqual(tc,full(sum(U,1)),ones(1,n));
        verifyEqual(tc,full(sum(U,2)),ones(n,1));
        p=zeros(n,1);
        for x=0:m-1
            for y=0:m-1
                j=1+x+m*y;
                v=mod(M*[x;y]+shift,m);
                p(j)=1+v(1)+m*v(2);
            end
        end
        verifyEqual(tc,d.permutation,p);
        verifyEqual(tc,full(U),full(sparse((1:n)',p,1,n,n)));
        verifyEqual(tc,d.inverse(p),(1:n)');
        verifyEqual(tc,sort(d.permutation),(1:n)');
        verifyEqual(tc,sum(d.cycle_lengths),n);
        verifyEqual(tc,numel(d.cycles),numel(d.cycle_lengths));
        verifyEqual(tc,sum(d.cycle_lengths==1),sum(p==(1:n)'));
    end
end
end

function test_exact_cycle_eigenphases_and_fixed_point_moments(tc)
for m=[2 3 4 5 6]
    M=[2 1;1 1];
    U=finite_koopman_torus(m,M);
    S=finite_koopman_torus(m,M,'spectral');
    D=finite_koopman_torus(m,M,'map');
    n=m*m;
    z=S.finite.eigenvalues();
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,numel(z),n);
    verifyLessThan(tc,max(abs(abs(z)-1)),1e-12);
    actual=eig(full(U));
    verifyLessThan(tc,max(arrayfun(@(x)min(abs(z-x)),actual)), ...
        2e-8);
    verifyEqual(tc,S.cycle_lengths,D.cycle_lengths);
    for k=[-7 -3 -1 0 1 2 3 5 9]
        fixedCount=0;
        for j=1:numel(D.cycles)
            ell=numel(D.cycles{j});
            if mod(k,ell)==0,fixedCount=fixedCount+ell;end
        end
        moment=S.finite.moments(k);
        verifyLessThan(tc,abs(moment-fixedCount/n),1e-14);
        verifyLessThan(tc,abs(moment-mean(z.^k)),6e-13);
        verifyLessThan(tc,abs(moment-trace(full(U)^k)/n), ...
            1e-10);
    end
end
end

function test_dual_fourier_character_permutation_and_translation(tc)
% The 2D character basis is built *independently*, avoiding the
% generator's dual-permutation formula.
for m=[2 3 4 5]
    n=m*m;
    x=mod((0:n-1)',m);y=floor((0:n-1)'/m);
    C=exp(2i*pi*(x*x'+y*y')/m)/sqrt(n);
    verifyLessThan(tc,norm(C'*C-eye(n),'fro'),2e-13*n);
    for model={struct('M',[2 1;1 1],'b',[0;0]), ...
               struct('M',[1 1;0 1],'b',[0;0]), ...
               struct('M',[0 -1;1 0],'b',[1;2])}
        M=model{1}.M;shift=model{1}.b;
        U=finite_koopman_torus(m,M,'matrix',shift);
        D=finite_koopman_torus(m,M,'map',shift);
        dual=C'*full(U)*C;
        V=sparse(D.dual_frequency_permutation,(1:n)', ...
            D.dual_frequency_phases,n,n);
        verifyLessThan(tc,norm(dual-full(V),'fro'),3e-12*n);
        verifyLessThan(tc,norm(full(U)*C-C*full(V),'fro'),3e-12*n);
    end
end
end

function test_operator_adjoint_multiple_RHS_and_size(tc)
for m=[2 5 11 31]
    n=m*m;
    M=[2 1;1 1];
    O=finite_koopman_torus(m,M,'operator',[3;-2]);
    verifyEqual(tc,O([],'size'),[n,n]);
    X=reshape(sin(.43*(1:3*n)),n,3) ...
      +1i*reshape(cos(.21*(1:3*n)),n,3);
    Y=O(X,'notransp');
    verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'),1e-14*n);
    verifyLessThan(tc,norm(O(X,'transp')-O(X,'ctransp'),'fro'), ...
        1e-14*n);
    verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')),1e-14*n);
    if m<=11
        U=finite_koopman_torus(m,M,'matrix',[3;-2]);
        verifyLessThan(tc,norm(Y-U*X,'fro'),1e-13*n);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            1e-13*n);
    end
    verifyError(tc,@()O(zeros(n-1,1),'notransp'), ...
        'um_finite_koopman_apply:Dimension');
    verifyError(tc,@()O(X,'other'), ...
        'um_finite_koopman_apply:Flag');
end
end

function test_translation_shift_cycles_and_permutation_order(tc)
for m=[2 3 4 7]
    D=finite_koopman_torus(m,eye(2),'map',[1;0]);
    verifyEqual(tc,sort(D.cycle_lengths),m*ones(m,1));
    verifyEqual(tc,D.order,m);
    U=finite_koopman_torus(m,eye(2));
    verifyEqual(tc,full(U),eye(m*m));
    D=finite_koopman_torus(m,eye(2),'map',[0;0]);
    verifyEqual(tc,D.cycle_lengths,ones(m*m,1));
end
end

function test_properties_inputs_and_nonbijective_rejection(tc)
[A,p]=finite_koopman_torus();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'orthogonal')));
verifyTrue(tc,any(strcmp(p,'sparse')));
verifyEqual(tc,size(finite_koopman_torus(3)),[9 9]);
verifyError(tc,@()finite_koopman_torus(1), ...
    'finite_koopman_torus:Order');
verifyError(tc,@()finite_koopman_torus(4,[2 0;0 1]), ...
    'finite_koopman_torus:Invertibility');
verifyEqual(tc,size(finite_koopman_torus(3,[2 0;0 1])),[9 9]);
verifyError(tc,@()finite_koopman_torus(4,[1 .5;0 1]), ...
    'finite_koopman_torus:Matrix');
verifyError(tc,@()finite_koopman_torus(4,[],[],'bad'), ...
    'finite_koopman_torus:Translation');
verifyError(tc,@()finite_koopman_torus(4,[],'oops'), ...
    'finite_koopman_torus:Mode');
verifyError(tc,@()finite_koopman_torus(65), ...
    'finite_koopman_torus:ExplicitOrder');
end
