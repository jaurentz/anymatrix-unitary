function tests=test_extended_lipkin_propagator
%TEST_EXTENDED_LIPKIN_PROPAGATOR Numerical and analytical model checks.
tests=functiontests(localfunctions);
end

function test_hamiltonian_independent_schwinger_construction(tc)
for n=[1 2 3 6 9]
    for alpha=[0 .2 -sqrt(.5)]
        for lambda=[0 .25 1]
            H=extended_lipkin_propagator(n,lambda,alpha,.3,'hamiltonian');
            k=(0:n)';
            S=n/2;
            sz=spdiags(k-S,0,n+1,n+1);
            sp=sparse(2:n+1,1:n,sqrt((1:n)'.*(n-(1:n)'+1)), ...
                n+1,n+1);
            sx=(sp+sp')/2;
            n_t=S*speye(n+1)+sz;
            Q=2*sx+alpha*n_t;
            expected=(1-lambda)*n_t-(lambda/n)*Q^2;
            verifyLessThan(tc,norm(full(H-expected),'fro'),2e-12*n);
            verifyLessThan(tc,norm(full(H-H'),'fro'),2e-12*n);
            [i,j]=find(H);
            verifyTrue(tc,all(abs(i-j)<=2));
            if alpha==0
                verifyLessThan(tc,norm(full(H(1:2:end,2:2:end)),'fro'), ...
                    2e-12*n);
            end
        end
    end
end
end

function test_evolution_unitarity(tc)
for n=[1 2 6 15]
    for lambda=[0 .2 .7 1]
        for alpha=[0 -.35 .6]
            U=extended_lipkin_propagator(n,lambda,alpha,.71);
            H=extended_lipkin_propagator(n,lambda,alpha,.71,'hamiltonian');
            verifyEqual(tc,size(U),[n+1,n+1]);
            verifyLessThan(tc,norm(U'*U-eye(n+1),'fro'),3e-12*n);
            verifyLessThan(tc,norm(U-expm(-.71i*full(H)),'fro'), ...
                2e-12*n);
        end
    end
end
end

function test_free_bosons_and_x_squared_exact_eigenphases(tc)
for n=[1 2 3 6 9]
    for t=[0 .43 1.6]
        A=extended_lipkin_propagator(n,0,.8,t);
        z=exp(-1i*t*(0:n)');
        verifyLessThan(tc,norm(A-diag(z),'fro'),2e-12*n);
        S=extended_lipkin_propagator(n,0,.8,t,'spectral');
        verifyEqual(tc,S.finite.kind,'exact');
        verifyLessThan(tc,abs(trace(A)-sum(S.finite.eigenvalues())), ...
            2e-12*n);
        B=extended_lipkin_propagator(n,1,0,t);
        T=extended_lipkin_propagator(n,1,0,t,'spectral');
        w=T.finite.eigenvalues();
        verifyEqual(tc,numel(w),n+1);
        verifyLessThan(tc,max(arrayfun(@(v)min(abs(v-w)),eig(B))), ...
            3e-11*n);
        verifyLessThan(tc,abs(trace(B)-sum(w)),5e-11*n);
        for p=[-2 0 1 3]
            verifyLessThan(tc,abs(trace(B^p)/(n+1)-T.finite.moments(p)), ...
                3e-11*n);
        end
    end
end
end

function test_critical_parameter_and_unknown_spectral_mode(tc)
n=16;alpha=sqrt(.5);lambda=1/(5+alpha^2);
S=extended_lipkin_propagator(n,lambda,alpha,1,'spectral');
verifyEqual(tc,S.family,'extended_lipkin_propagator');
verifyEqual(tc,S.dimension,n+1);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEqual(tc,S.limit.kind,'unknown');
U=extended_lipkin_propagator(n,lambda,alpha);
verifyLessThan(tc,norm(U'*U-eye(n+1),'fro'),1e-11);
A=extended_lipkin_propagator(n,lambda,alpha,0,'spectral');
verifyEqual(tc,A.finite.kind,'exact');
verifyEqual(tc,A.finite.eigenvalues(),ones(n+1,1));
end

function test_anymatrix_contract_validation(tc)
[e,p]=extended_lipkin_propagator();
verifyEmpty(tc,e);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'square')));
verifyError(tc,@()extended_lipkin_propagator(0), ...
    'MATLAB:expectedPositive');
verifyError(tc,@()extended_lipkin_propagator(5,1.1), ...
    'MATLAB:expectedLessEqual');
verifyError(tc,@()extended_lipkin_propagator(4,.2,.4,1,'operator'), ...
    'extended_lipkin_propagator:Mode');
verifyError(tc,@()extended_lipkin_propagator(513,.2,.4,1), ...
    'extended_lipkin_propagator:ExplicitOrder');
end
