function tests=test_floquet_extended
%TEST_FLOQUET_EXTENDED Tests for rotor, dual-unitary and XXZ Floquet models.
tests=functiontests(localfunctions);
end

function test_kicked_rotor_analytic_limits(testCase)
for n=[2 3 6 9]
    for hbar=[0.43 2*pi/n]
        q=(0:n-1)';
        D=diag(exp(-1i*hbar*q.^2/2));
        U=kicked_rotor(n,0,hbar);
        verifyEqual(testCase,size(U),[n n]);
        verifyLessThan(testCase,norm(U-D,'fro'),2e-12*n);
        check_unitarity(testCase,U);
        for shift=[0 .27 .8]
            V=kicked_rotor(n,4.7,hbar,shift);
            check_unitarity(testCase,V);
            verifyLessThan(testCase,norm(V-kicked_rotor(n,4.7,hbar,shift),'fro'),1e-15);
        end
    end
end
% Explicit two-dimensional rotor with Hadamard Fourier matrix:
% F'*diag(e^{-i A},e^{i A})*F = [cos A,-i sin A;-i sin A,cos A].
hbar=.75;K=1.3;A=K/hbar;
expected=diag([1;exp(-1i*hbar/2)])* ...
         [cos(A),-1i*sin(A);-1i*sin(A),cos(A)];
verifyLessThan(testCase,norm(kicked_rotor(2,K,hbar,0)-expected,'fro'),3e-14);
thrown=false;
try,kicked_rotor(2,1,0);catch,thrown=true;end
verifyTrue(testCase,thrown);
end

function test_two_site_embedding_and_order(testCase)
Swap=[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
for n=[2 3 4]
    for a=1:n
        for b=1:n
            if a==b,continue;end
            U=qm_two_site_gate(n,Swap,a,b);
            N=2^n;
            reference=speye(N);
            reference=reference*0;
            for k=0:N-1
                s=uint64(k);
                x=bitget(s,n-a+1);
                y=bitget(s,n-b+1);
                out=bitset(bitset(s,n-a+1,y),n-b+1,x);
                reference(double(out)+1,k+1)=1;
            end
            verifyLessThan(testCase,norm(U-reference,'fro'),1e-13);
            check_unitarity(testCase,U);
        end
    end
end
end

function test_dual_unitarity_of_local_gate(testCase)
for chi=[0 .31 pi/8 .71]
    for angles={[0 0 0 0],[.37 -.61 .29 .83]}
        G=dual_unitary_circuit(2,1,chi,angles{1},'open');
        check_unitarity(testCase,G);
        T=zeros(4);
        for oa=0:1
            for ob=0:1
                for ia=0:1
                    for ib=0:1
                        T(2*oa+ia+1,2*ob+ib+1)=G(2*oa+ob+1,2*ia+ib+1);
                    end
                end
            end
        end
        check_unitarity(testCase,T);
    end
end
Swap=[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
verifyLessThan(testCase,norm(dual_unitary_circuit(2,1,0,[0 0 0 0],'open')-Swap,'fro'),1e-14);
end

function test_dual_circuit_brickwork(testCase)
for n=[2 3 4 6]
    if mod(n,2),boundary='open';else,boundary='periodic';end
    if n==2,boundary='open';end
    U=dual_unitary_circuit(n,1,.31,[.2 -.3 .4 -.5],boundary);
    check_unitarity(testCase,U);
    check_unitarity(testCase,dual_unitary_circuit(n));
    verifyLessThan(testCase,norm(dual_unitary_circuit(n,2,.31,[.2 -.3 .4 -.5],boundary)-U^2,'fro'),2e-12*2^n);
    I=dual_unitary_circuit(n,0,.31,[.2 -.3 .4 -.5],boundary);
    verifyLessThan(testCase,norm(I-eye(2^n),'fro'),1e-14);
end
verifyError(testCase,@()dual_unitary_circuit(5,1,.3,[0 0 0 0],'periodic'), ...
    'dual_unitary_circuit:Boundary');
end

function test_xxz_unitary_r_and_yang_baxter(testCase)
for gamma=[.3 .8 1.4]
    for u=[-.5 0 .31 .8]
        G=xxz_r_gate(u,gamma);
        check_unitarity(testCase,G);
        if u==0
            verifyLessThan(testCase,norm(G-1i*eye(4),'fro'),2e-14);
        end
        % The gate preserves excitation number in its 2-qubit basis.
        verifyLessThan(testCase,norm(G([1,4],[2,3]),'fro')+ ...
                              norm(G([2,3],[1,4]),'fro'),1e-14);
    end
    u=.28;v=-.19;
    L=kron(xxz_r_gate(u,gamma),eye(2))* ...
      kron(eye(2),xxz_r_gate(u+v,gamma))* ...
      kron(xxz_r_gate(v,gamma),eye(2));
    R=kron(eye(2),xxz_r_gate(v,gamma))* ...
      kron(xxz_r_gate(u+v,gamma),eye(2))* ...
      kron(eye(2),xxz_r_gate(u,gamma));
    verifyLessThan(testCase,norm(L-R,'fro'),2e-12);
end
end

function test_xxz_floquet_sectors_and_scalar_limit(testCase)
for n=[2 3 4 5 6]
    if mod(n,2) || n==2,boundary='open';else,boundary='periodic';end
    U=xxz_floquet(n,.42,.8,boundary);
    check_unitarity(testCase,U);
    check_unitarity(testCase,xxz_floquet(n));
    N=2^n;
    state=uint64((0:N-1)');
    excit=zeros(N,1);
    for k=1:n,excit=excit+double(bitget(state,k));end
    [row,col]=find(abs(U)>1e-10);
    verifyTrue(testCase,all(excit(row)==excit(col)));
    numgates=n-1;
    if strcmp(boundary,'periodic'),numgates=n;end
    phase=1i^numgates;
    V=xxz_floquet(n,0,.8,boundary);
    verifyLessThan(testCase,norm(V-phase*eye(N),'fro'),2e-11*N);
end
verifyError(testCase,@()xxz_floquet(5,.4,.8,'periodic'),'xxz_floquet:Boundary');
thrown=false;
try,xxz_floquet(4,.4,0,'periodic');catch,thrown=true;end
verifyTrue(testCase,thrown);
end

function check_unitarity(testCase,U)
n=size(U,1);
verifyEqual(testCase,size(U),[n,n]);
verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),8e-11*n);
end
