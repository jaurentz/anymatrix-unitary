function tests=test_butson_hadamard
tests=functiontests(localfunctions);
end

function test_tao_exponents_and_exact_orthogonality(tc)
R=butson_hadamard(6,'tao','exponents');
verifyEqual(tc,R.q,3);
verifyEqual(tc,size(R.E),[6 6]);
verifyEqual(tc,unique(R.E),[0;1;2]);
for i=1:6
    for j=1:i-1
        delta=mod(R.E(i,:)-R.E(j,:),3);
        verifyEqual(tc,[sum(delta==0),sum(delta==1),sum(delta==2)], ...
            [2 2 2]);
    end
end
U=butson_hadamard(6,'tao');
verifyLessThan(tc,norm(U'*U-eye(6),'fro'),3e-14);
verifyLessThan(tc,norm(abs(U)-ones(6)/sqrt(6),'fro'),3e-14);
verifyLessThan(tc,norm(U-U.','fro'),3e-14);
S=butson_hadamard(6,'tao','spectral');
verifyEqual(tc,S.finite.kind,'unknown');
% Left Hadamard equivalence preserves the Gram matrix but not eigenphases.
V=diag(exp(1i*[.3;.2;-.4;.1;.5;-.1]))*U;
verifyLessThan(tc,norm(V'*V-eye(6),'fro'),3e-14);
verifyGreaterThan(tc,abs(trace(V)-trace(U)),.05);
end

function test_paley_orders_and_legendre_cores(tc)
for p=[3 7 11 19 23]
    n=p+1;
    R=butson_hadamard(n,'paley','exponents');
    verifyEqual(tc,R.q,2);
    H=1-2*R.E;
    verifyTrue(tc,all(abs(H(:))==1));
    verifyEqual(tc,H*H',n*eye(n)); % exact integer cancellation
    verifyEqual(tc,H(1,:),ones(1,n));
    verifyEqual(tc,H(:,1),ones(n,1));
    Q=H(2:end,2:end)+eye(p);
    verifyEqual(tc,Q+Q',zeros(p)); % Legendre oddness
    U=butson_hadamard(n,'paley');
    verifyTrue(tc,isreal(U));
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),5e-13*n);
    verifyLessThan(tc,norm(abs(U)-ones(n)/sqrt(n),'fro'),5e-13*n);
    S=butson_hadamard(n,'paley','spectral');
    verifyEqual(tc,S.finite.kind,'unknown');
end
end

function test_independent_qth_roots_and_signs(tc)
T=butson_hadamard(6,'tao','exponents');
U=butson_hadamard(6,'tao');
verifyLessThan(tc,norm(sqrt(6)*U-exp(2i*pi*T.E/3),'fro'),1e-14);
P=butson_hadamard(12,'paley','exponents');
V=butson_hadamard(12,'paley');
verifyEqual(tc,sqrt(12)*V,1-2*P.E,'AbsTol',2e-14);
end

function test_unsupported_orders_and_modes(tc)
verifyError(tc,@()butson_hadamard(12,'tao'),'butson_hadamard:Order');
verifyError(tc,@()butson_hadamard(10,'paley'),'butson_hadamard:Order');
verifyError(tc,@()butson_hadamard(13,'paley'),'butson_hadamard:Order');
verifyError(tc,@()butson_hadamard(6,'bogus'), ...
    'butson_hadamard:Construction');
verifyError(tc,@()butson_hadamard(6,'tao','bad'), ...
    'butson_hadamard:Mode');
[U,props]=butson_hadamard();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(props,'unitary')));
end
