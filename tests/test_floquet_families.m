function tests=test_floquet_families
%TEST_FLOQUET_FAMILIES Independent unitary and spectral regression checks.
tests=functiontests(localfunctions);
end

function test_kicked_top_limits(testCase)
for n=[2 3 4 7]
    j=(n-1)/2;
    m=(j:-1:-j)';
    for kappa=[0 0.8 3.2]
        U=kicked_top(n,kappa,0);
        check_unitary(testCase,U);
        expected=exp(-1i*kappa*m.^2/(2*j));
        verifyLessThan(testCase,norm(U-diag(expected),'fro'),3e-13*n);
    end
    p=0.37;
    U=kicked_top(n,0,p);
    check_unitary(testCase,U);
    for k=1:4
        verifyLessThan(testCase,abs(trace(U^k)-sum(exp(-1i*k*p*m))),2e-11*n);
    end
    check_unitary(testCase,kicked_top(n,3.1,.73));
end
end

function test_floquet_ising_exact_layers(testCase)
for n=[2 3 4 5]
    N=2^n;
    for boundary={'open','periodic'}
        U=floquet_ising(n,.43,.62,boundary{1});
        verifyEqual(testCase,size(U),[N N]);
        check_unitary(testCase,U);
        check_unitary(testCase,floquet_ising(n,pi/4,pi/4,boundary{1}));
        verifyLessThan(testCase,norm(floquet_ising(n,0,0,boundary{1})-eye(N),'fro'),2e-12);
        D=floquet_ising(n,.43,0,boundary{1});
        verifyLessThan(testCase,norm(D-diag(diag(D)),'fro'),1e-12);
        bits=uint64((0:N-1)');
        z=zeros(N,n);
        for j=1:n
            z(:,j)=1-2*double(bitget(bits,j));
        end
        energy=sum(z(:,1:end-1).*z(:,2:end),2);
        if strcmp(boundary{1},'periodic')
            energy=energy+z(:,end).*z(:,1);
        end
        verifyLessThan(testCase,norm(diag(D)-exp(-1i*.43*energy)),2e-12);
    end
    R=[cos(.33),-1i*sin(.33);-1i*sin(.33),cos(.33)];
    X=1;
    for j=1:n,X=kron(X,R);end
    verifyLessThan(testCase,norm(floquet_ising(n,0,.33)-X,'fro'),3e-12);
end
verifyError(testCase,@()floquet_ising(3,.4,.5,'bad'),'floquet_ising:Boundary');
end

function test_selfdual_and_integrable_comparison(testCase)
for n=[3 4 5]
    U=selfdual_kicked_ising(n,zeros(n,1));
    V=floquet_ising(n,pi/4,pi/4,'periodic');
    check_unitary(testCase,U);
    for k=1:5
        verifyLessThan(testCase,abs(trace(U^k)-trace(V^k)),4e-10*2^n);
    end
    U=selfdual_kicked_ising(n,0.37);
    V=selfdual_kicked_ising(n,0.37*ones(1,n));
    verifyEqual(testCase,U,V,'AbsTol',2e-13);
    check_unitary(testCase,selfdual_kicked_ising(n,0.2*sin(1:n)));
end
verifyError(testCase,@()selfdual_kicked_ising(4,ones(1,3)), ...
    'selfdual_kicked_ising:Fields');
end

function test_driven_kitaev_ph_symmetry(testCase)
for n=[2 3 5 7]
    for boundary={'open','periodic'}
        if n==2 && strcmp(boundary{1},'periodic'),continue;end
        opts=struct('mu1',.2,'mu2',2.4,'hopping',1.1, ...
                    'pairing',.9,'tau1',.4,'tau2',.7, ...
                    'boundary',boundary{1});
        U=driven_kitaev(n,opts);
        verifyEqual(testCase,size(U),[2*n 2*n]);
        check_unitary(testCase,U);
        C=[zeros(n),eye(n);eye(n),zeros(n)];
        verifyLessThan(testCase,norm(C*conj(U)*C-U,'fro'),5e-11*n);
    end
    verifyLessThan(testCase,norm(driven_kitaev(n,struct('tau1',0,'tau2',0)) ...
          -eye(2*n),'fro'),2e-13);
end
% Exactly soluble static open chain: |hopping|=|pairing|, mu=0,
% produces two zero-energy Majorana modes (two eigenvalues at +1).
opts=struct('mu1',0,'mu2',0,'hopping',1,'pairing',1, ...
            'tau1',.23,'tau2',.51,'boundary','open');
U=driven_kitaev(6,opts);
verifyGreaterThanOrEqual(testCase,sum(abs(eig(U)-1)<1e-10),2);
verifyError(testCase,@()driven_kitaev(3,struct('nonsense',4)), ...
    'driven_kitaev:Options');
verifyError(testCase,@()driven_kitaev(2,struct('boundary','periodic')), ...
    'driven_kitaev:Boundary');
end

function check_unitary(testCase,U)
n=size(U,1);
verifyEqual(testCase,size(U),[n,n]);
verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),5e-11*n);
end
