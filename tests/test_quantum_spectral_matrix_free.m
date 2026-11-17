function tests=test_quantum_spectral_matrix_free
%TEST_QUANTUM_SPECTRAL_MATRIX_FREE Analytic laws and large tensor actions.
tests=functiontests(localfunctions);
end

function test_unimodular_densities_and_moments(testCase)
names={'haar_special_unitary','haar_unimodular_coe','haar_unimodular_cse'};
for j=1:numel(names)
    f=names{j};
    for n=2:7
        S=feval(f,n,[],'spectral');
        verifyEqual(testCase,S.ensemble.kind,'exact_expectation');
        verifyEqual(testCase,S.dimension,n*(1+strcmp(f,'haar_unimodular_cse')));
        rho=S.ensemble.density;
        verifyLessThan(testCase,abs(integral(rho,-pi,pi)-1),3e-11);
        grid=linspace(-pi,pi,97);
        verifyGreaterThanOrEqual(testCase,min(rho(grid)),-1e-12);
        verifyLessThan(testCase,norm(rho(grid)-rho(grid+2*pi)),3e-11);
        verifyEqual(testCase,S.ensemble.moments(0),1,'AbsTol',1e-14);
        for k=[1 2 3 n n+1 2*n 3*n]
            m=S.ensemble.moments(k);
            % Test Fourier moments against independent phase quadrature.
            theta=linspace(-pi,pi,4097);
            theta=theta(1:end-1);
            if strcmp(f,'haar_unimodular_coe') && n==2
                integrand=@(q)rho(q).*exp(1i*k*q);
                mk=integral(integrand,-pi,0)+integral(integrand,0,pi);
            else
                mk=2*pi/numel(theta)*sum(rho(theta).*exp(1i*k*theta));
            end
            verifyLessThan(testCase,abs(m-mk),4e-10);
            verifyLessThan(testCase, ...
                abs(S.ensemble.moments(-k)-conj(m)),2e-13);
        end
    end
end
end

function test_coe_exception_and_order_one_atoms(testCase)
S=haar_unimodular_coe(2,[],'spectral');
theta=linspace(-pi,pi,81);
verifyEqual(testCase,S.ensemble.density(theta),abs(sin(theta))/4, ...
    'AbsTol',1e-14);
verifyLessThan(testCase,abs(S.ensemble.moments(2)+1/3),1e-13);
verifyLessThan(testCase,abs(S.ensemble.moments(4)+1/15),1e-13);
verifyEqual(testCase,S.ensemble.moments(3),0);
for f={@haar_special_unitary,@haar_unimodular_coe,@haar_unimodular_cse}
    T=f{1}(1,[],'spectral');
    verifyEqual(testCase,T.ensemble.kind,'atomic');
    verifyEqual(testCase,T.finite.kind,'exact');
    verifyEqual(testCase,unique(T.finite.eigenvalues()),1);
end
end

function test_entangled_reference_exact_without_eig(testCase)
for d=1:5
    for phi=[0 pi .37]
        S=entangled_reflection(d,phi,'spectral');
        A=entangled_reflection(d,phi);
        verifyEqual(testCase,S.dimension,d^2);
        verifyEqual(testCase,S.finite.kind,'exact');
        for k=-3:4
            verifyLessThan(testCase,abs(S.finite.moments(k) ...
                 -trace(A^k)/d^2),1e-12);
        end
        z=S.finite.eigenvalues();
        a=eig(full(A));
        verifyEqual(testCase,numel(z),d^2);
        verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),a)),3e-12);
    end
end
end

function test_tensor_permutation_spectral_and_operator(testCase)
for d=[2 3]
    for perm={[2 1],[2 3 1],[3 1 2]}
        p=perm{1};
        P=tensor_permutation(d,p);
        S=tensor_permutation(d,p,'spectral');
        O=tensor_permutation(d,p,'operator');
        verifyEqual(testCase,O([],'size'),size(P));
        verifyEqual(testCase,S.dimension,size(P,1));
        X=randn(size(P,1),3)+1i*randn(size(P,1),3);
        verifyLessThan(testCase,norm(O(X,'notransp')-P*X,'fro'),3e-12);
        verifyLessThan(testCase,norm(O(X,'ctransp')-P'*X,'fro'),3e-12);
        verifyLessThan(testCase,norm(O(X,'transp')-P'*X,'fro'),3e-12);
        for k=-4:4
            verifyLessThan(testCase,abs(S.finite.moments(k) ...
                -trace(P^k)/size(P,1)),3e-12);
        end
        z=S.finite.eigenvalues();
        verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),eig(full(P)))),3e-12);
    end
end
% O(2^12) applies a tensor reshuffle with no sparse matrix stored.
P=tensor_permutation(2,[12 2:11 1],'operator');
verifyEqual(testCase,P([],'size'),[4096 4096]);
x=randn(4096,1)+1i*randn(4096,1);
verifyLessThan(testCase,norm(P(P(x,'notransp'),'ctransp')-x),1e-12);
end

function test_haar_tensor_power_spectral_and_operator(testCase)
state=rng;
d=2;t=4;seed=31;
A=haar_tensor_power(d,t,seed);
O=haar_tensor_power(d,t,seed,'operator');
S=haar_tensor_power(d,t,seed,'spectral');
verifyEqual(testCase,O([],'size'),[16 16]);
checkactions(testCase,A,O);
verifyEqual(testCase,S.finite.kind,'exact');
for k=-3:4
    verifyLessThan(testCase,abs(S.finite.moments(k)-trace(A^k)/16),2e-12);
end
z=S.finite.eigenvalues();
verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),eig(A))),3e-12);
verifyEqual(testCase,rng,state);
verifyError(testCase,@()haar_tensor_power(2,3,[],'spectral'), ...
    'unitary_spectral_reference:Seed');
% The prior dense cap is intentionally bypassed in operator mode.
H=haar_tensor_power(2,13,seed,'operator');
verifyEqual(testCase,H([],'size'),[8192 8192]);
x=randn(8192,2)+1i*randn(8192,2);
verifyLessThan(testCase,norm(H(H(x,'notransp'),'ctransp')-x,'fro'),2e-10);
end

function test_independent_block_tensor_reference(testCase)
state=rng;
A=haar_block_tensor_power(2,2,1,47);
O=haar_block_tensor_power(2,2,1,47,'operator');
S=haar_block_tensor_power(2,2,1,47,'spectral');
checkactions(testCase,A,O);
verifyEqual(testCase,S.dimension,8);
for k=-3:4
    verifyLessThan(testCase,abs(S.finite.moments(k)-trace(A^k)/8),2e-12);
end
z=S.finite.eigenvalues();
verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),eig(A))),3e-12);
verifyEqual(testCase,rng,state);
verifyError(testCase,@()haar_block_tensor_power(2,2,1,[],'spectral'), ...
    'unitary_spectral_reference:Seed');
H=haar_block_tensor_power(2,7,6,47,'operator');
verifyEqual(testCase,H([],'size'),[8192 8192]);
x=randn(8192,1)+1i*randn(8192,1);
verifyLessThan(testCase,norm(H(H(x,'notransp'),'ctransp')-x),2e-10);
end

function checkactions(testCase,A,O)
n=size(A,1);
X=randn(n,3)+1i*randn(n,3);
verifyLessThan(testCase,norm(O(X,'notransp')-A*X,'fro'),2e-11*n);
verifyLessThan(testCase,norm(O(X,'transp')-A'*X,'fro'),2e-11*n);
end
