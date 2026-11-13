function tests=test_quantum_expansion
%TEST_QUANTUM_EXPANSION Algebraic and spectral checks of three physical models.
tests=functiontests(localfunctions);
end

function test_bose_hubbard_exact_free_hopping(testCase)
for n=[0 1 2 3 7 11]
    for J=[0 .6 1.1]
        t=.71;
        U=kicked_bose_hubbard(n,J,0,t,1,0);
        check_unitary(testCase,U);
        verifyEqual(testCase,size(U),[n+1,n+1]);
        energies=J*(2*(0:n)'-n);
        for k=1:4
            verifyLessThan(testCase, ...
                abs(trace(U^k)-sum(exp(-1i*k*t*energies))),3e-11*(n+1));
        end
    end
end
end

function test_bose_hubbard_diagonal_kick_and_expm(testCase)
for n=[0 1 2 4 8]
    m=(0:n)';
    J=.65;g=.43;bias=.27;t1=.62;t2=.81;
    energy=.5*g*(m.*(m-1)+(n-m).*(n-m-1)) + .5*bias*(2*m-n);
    D=diag(exp(-1i*t2*energy));
    V=kicked_bose_hubbard(n,0,g,t1,t2,bias);
    verifyLessThan(testCase,norm(V-D,'fro'),1e-12*(n+1));
    W=kicked_bose_hubbard(n,J,g,t1,t2,bias);
    check_unitary(testCase,W);
    if n>0
        hop=sqrt(((1:n).*(n-(1:n)+1))).';
        H=-J*(diag(hop,1)+diag(hop,-1));
    else
        H=0;
    end
    reference=D*expm(-1i*t1*H);
    verifyLessThan(testCase,norm(W-reference,'fro'),2e-12*(n+1));
    verifyLessThan(testCase,norm(W-kicked_bose_hubbard(n,J,g,t1,t2,bias),'fro'),1e-14);
end
verifyLessThan(testCase,norm(kicked_bose_hubbard(3,1,.5,0,0)-eye(4),'fro'),1e-14);
verifyTrue(testCase,raises(@()kicked_bose_hubbard(-1)));
end

function test_jaynes_cummings_exact_spectra(testCase)
for n=[0 1 2 4 9]
    for g=[0 .8]
        for delta=[0 .37]
            t=.65;wc=1.2;wa=wc+delta;
            U=jaynes_cummings(n,g,delta,t,wc);
            check_unitary(testCase,U);
            verifyTrue(testCase,issparse(U));
            verifyEqual(testCase,size(U),[2*(n+1),2*(n+1)]);
            expected=[exp(1i*t*wa/2);exp(-1i*t*(wc*n+wa/2))];
            for k=0:n-1
                center=wc*(k+.5);
                split=hypot(delta/2,g*sqrt(k+1));
                expected=[expected;exp(-1i*t*(center+split)); ...
                                    exp(-1i*t*(center-split))]; %#ok<AGROW>
            end
            for r=1:5
                verifyLessThan(testCase,abs(trace(U^r)-sum(expected.^r)), ...
                   3e-11*(n+1));
            end
        end
    end
end
end

function test_jaynes_cummings_matches_truncated_H(testCase)
for n=[0 1 2 5]
    for g=[0 .51]
        delta=-.27;wc=1.3;t=.43;
        wa=wc+delta;
        d=2*(n+1);
        H=zeros(d);
        for k=0:n
            H(2*k+1,2*k+1)=wc*k-wa/2;
            H(2*k+2,2*k+2)=wc*k+wa/2;
        end
        for k=0:n-1
            i=2*k+2;j=i+1;
            H(i,j)=g*sqrt(k+1);
            H(j,i)=H(i,j);
        end
        U=jaynes_cummings(n,g,delta,t,wc);
        reference=expm(-1i*t*H);
        verifyLessThan(testCase,norm(full(U)-reference,'fro'),4e-12*d);
    end
end
verifyLessThan(testCase,norm(full(jaynes_cummings(4,0,0,0))-eye(10),'fro'),1e-13);
verifyLessThan(testCase,norm(full(jaynes_cummings(4,0,0,1,0))-eye(10),'fro'),1e-13);
verifyTrue(testCase,raises(@()jaynes_cummings(2,1,0,Inf)));
end

function test_chalker_coddington_sparse_unitary_and_pi_pairing(testCase)
for grid=[1 1; 2 1; 1 3; 2 3; 4 2]'
    nx=grid(1);ny=grid(2);
    N=4*nx*ny;
    phases=.23*(1:N)'+.4*sin((1:N)');
    for theta={[0 0], [.19 1.03], [pi/4 pi/4], [pi/2 pi/2]}
        U=chalker_coddington(nx,ny,theta{1},phases);
        check_unitary(testCase,U);
        verifyTrue(testCase,issparse(U));
        verifyEqual(testCase,size(U),[N N]);
        verifyLessThanOrEqual(testCase,nnz(U),2*N);
        x=mod((0:N-1)',2*nx);
        y=floor((0:N-1)'/(2*nx));
        C=spdiags((-1).^(x+y),0,N,N);
        verifyLessThan(testCase,norm(C*U*C+U,'fro'),1e-13*N);
        for k=[1 3 5]
            verifyLessThan(testCase,abs(trace(U^k)),2e-10*N);
        end
    end
end
end

function test_chalker_coddington_exact_network_and_phases(testCase)
% Reference mapping for the smallest 2x2 periodic directed-link grid.
U=chalker_coddington(1,1,0,0);
expected=sparse([2 3 4 1],[1 4 2 3],[1 1 1 -1],4,4);
verifyEqual(testCase,U,expected);
verifyLessThan(testCase,norm(full(U^4+speye(4)),'fro'),1e-13);
for theta=[0 pi/2]
    U=chalker_coddington(2,3,theta);
    verifyTrue(testCase,all(full(sum(abs(U)>1e-12,1))==1));
    verifyTrue(testCase,all(full(sum(abs(U)>1e-12,2))==1));
    check_unitary(testCase,U);
end
theta=.4;phi=.27;
U=chalker_coddington(2,2,theta,0);
V=chalker_coddington(2,2,theta,phi);
verifyLessThan(testCase,norm(V-exp(1i*phi)*U,'fro'),1e-13);
phases=linspace(-pi,pi,16)';
Q=chalker_coddington(2,2,theta,phases);
D=spdiags(exp(1i*phases),0,16,16);
verifyLessThan(testCase,norm(Q-U*D,'fro'),1e-13);
for theta=[.31 pi/4]
    U=chalker_coddington(1,2,theta);
    verifyTrue(testCase,all(full(sum(spones(U),1))==2));
    verifyTrue(testCase,all(full(sum(spones(U),2))==2));
end
verifyError(testCase,@()chalker_coddington(2,3,[.1 .2 .3]), ...
    'chalker_coddington:Angles');
verifyError(testCase,@()chalker_coddington(2,3,pi/4,ones(3,1)), ...
    'chalker_coddington:Phases');
end

function check_unitary(testCase,U)
n=size(U,1);
verifyLessThan(testCase,norm(U'*U-speye(n),'fro'),6e-12*max(1,n));
end

function yes=raises(f)
yes=false;
try
    f();
catch
    yes=true;
end
end
