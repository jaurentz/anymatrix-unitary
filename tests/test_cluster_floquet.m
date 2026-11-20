function tests=test_cluster_floquet
%TEST_CLUSTER_FLOQUET Cluster symmetries, pi edge modes, exact spectra.
tests=functiontests(localfunctions);
end

function test_independent_cluster_hamiltonian_and_local_pulse(tc)
for n=[4 5 6]
    theta=.18+.12*(1:n)';
    variants={'open'};
    if mod(n,2)==0,variants{end+1}='periodic';end
    for boundary=variants
        for kick=[0 .19]
            A=cluster_floquet(n,theta,boundary{1},kick);
            N=2^n;
            verifyEqual(tc,size(A),[N N]);
            verifyLessThan(tc,norm(A'*A-eye(N),'fro'),5e-11*N);
            P=1;
            pulse=[-1i*sin(kick),cos(kick);cos(kick),-1i*sin(kick)];
            for j=1:n,P=kron(P,pulse);end
            H=zeros(N);
            if strcmp(boundary{1},'open'),sites=2:n-1;else,sites=1:n;end
            for j=sites
                operators=repmat({'I'},1,n);
                operators{j}='X';
                operators{mod(j-2,n)+1}='Z';
                operators{mod(j,n)+1}='Z';
                H=H+theta(j)*pauli_word(operators);
            end
            ref=P*expm(-1i*H);
            verifyLessThan(tc,norm(A-ref,'fro'),7e-11*N);
        end
    end
end
end

function test_Z2_x_Z2_symmetry_and_open_edge_anticonmutation(tc)
for n=[4 5 6]
    theta=.22+.08*sin(1:n);
    [Godd,Geven]=symmetries(n);
    for kick=[0 .18]
        U=cluster_floquet(n,theta,'open',kick);
        N=2^n;
        verifyLessThan(tc,norm(Godd*U-U*Godd,'fro'),4e-11*N);
        verifyLessThan(tc,norm(Geven*U-U*Geven,'fro'),4e-11*N);
        if kick==0
            left=pauli_word([{'Z'},repmat({'I'},1,n-1)]);
            right=pauli_word([repmat({'I'},1,n-1),{'Z'}]);
            verifyLessThan(tc,norm(U*left+left*U,'fro'),4e-11*N);
            verifyLessThan(tc,norm(U*right+right*U,'fro'),4e-11*N);
            otherLeft=pauli_word([{'X','Z'},repmat({'I'},1,n-2)]);
            otherRight=pauli_word([repmat({'I'},1,n-2),{'Z','X'}]);
            verifyLessThan(tc,norm(U*otherLeft+otherLeft*U,'fro'),4e-11*N);
            verifyLessThan(tc,norm(U*otherRight+otherRight*U,'fro'),4e-11*N);
        end
    end
end
for n=[4 6]
    U=cluster_floquet(n,.28+.03*(1:n),'periodic',0);
    [Godd,Geven]=symmetries(n);
    N=2^n;
    verifyLessThan(tc,norm(Godd*U-U*Godd,'fro'),4e-11*N);
    verifyLessThan(tc,norm(Geven*U-U*Geven,'fro'),4e-11*N);
    edge=pauli_word([{'Z'},repmat({'I'},1,n-1)]);
    verifyGreaterThan(tc,norm(U*edge+edge*U,'fro'),1e-2);
end
end

function test_paired_symmetry_sectors_at_pi(tc)
n=4;theta=[.17 .43 .89 .21];
U=cluster_floquet(n,theta,'open');
[Go,Ge]=symmetries(n);
Zleft=pauli_word({'Z','I','I','I'});
N=2^n;
for e=[-1 1]
    for o=[-1 1]
        proj=(eye(N)+e*Ge)*(eye(N)+o*Go)/4;
        partner=(eye(N)+e*Ge)*(eye(N)-o*Go)/4;
        verifyLessThan(tc,norm(U*proj-proj*U,'fro'),1e-12*N);
        verifyLessThan(tc,norm(partner*Zleft*proj-Zleft*proj,'fro'),1e-12*N);
        verifyLessThan(tc,norm(U*Zleft*proj+Zleft*U*proj,'fro'),1e-12*N);
        verifyEqual(tc,round(real(trace(proj))),4);
    end
end
end

function test_exact_finite_trace_spectra(tc)
for n=[4 5 6]
    theta=.13+.117*(1:n)';
    variants={'open'};
    if mod(n,2)==0,variants{end+1}='periodic';end
    for boundary=variants
        U=cluster_floquet(n,theta,boundary{1});
        S=cluster_floquet(n,theta,boundary{1},0,'spectral');
        verifyEqual(tc,S.dimension,2^n);
        verifyEqual(tc,S.finite.kind,'exact');
        z=S.finite.eigenvalues();
        verifyEqual(tc,numel(z),2^n);
        lam=eig(U);
        verifyLessThan(tc,max(arrayfun(@(w)min(abs(z-w)),lam)),5e-9);
        verifyLessThan(tc,max(arrayfun(@(w)min(abs(lam-w)),z)),5e-9);
        for k=-5:6
            verifyLessThan(tc,abs(mean(lam.^k)-S.finite.moments(k)), ...
                5e-10);
        end
        if strcmp(boundary{1},'open')
            for k=[-5 -3 -1 1 3 5]
                verifyEqual(tc,S.finite.moments(k),0,'AbsTol',1e-13);
            end
        end
    end
end
end

function test_bulk_disorder_vs_symmetry_preserving_kick(tc)
n=4;theta=[.21 .32 .43 .54];
ideal=cluster_floquet(n,theta,'open');
disordered=cluster_floquet(n,[3.1 -.82 .37 1.4],'open');
Z1=pauli_word({'Z','I','I','I'});
verifyLessThan(tc,norm(ideal*Z1+Z1*ideal,'fro'),1e-11);
verifyLessThan(tc,norm(disordered*Z1+Z1*disordered,'fro'),1e-11);
pert=cluster_floquet(n,theta,'open',.17);
[Godd,Geven]=symmetries(n);
verifyLessThan(tc,norm(pert*Godd-Godd*pert,'fro'),1e-11);
verifyLessThan(tc,norm(pert*Geven-Geven*pert,'fro'),1e-11);
verifyGreaterThan(tc,norm(pert*Z1+Z1*pert,'fro'),1e-2);
lambda=eig(pert);
deviation=max(arrayfun(@(z)min(abs(lambda+z)),lambda));
verifyGreaterThan(tc,deviation,1e-3); % pi-pairing no longer guaranteed
S=cluster_floquet(n,theta,'open',.17,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
end

function test_matrix_free_adjoint_and_validation(tc)
for n=[4 5 7]
    theta=.2+.05*(1:n);
    O=cluster_floquet(n,theta,'open',.21,'operator');
    N=2^n;
    A=cluster_floquet(n,theta,'open',.21);
    X=reshape(sin(.21*(1:3*N)),N,3)+ ...
        1i*reshape(cos(.31*(1:3*N)),N,3);
    verifyEqual(tc,O([],'size'),[N N]);
    verifyLessThan(tc,norm(O(X,'notransp')-A*X,'fro'),1e-10*N);
    verifyLessThan(tc,norm(O(X,'transp')-A'*X,'fro'),1e-10*N);
    verifyLessThan(tc,norm(O(X,'ctransp')-A'*X,'fro'),1e-10*N);
    verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
        1e-10*N);
end
n=14;N=2^n;
O=cluster_floquet(n,.35,'open',0,'operator');
X=complex(sin(.15*(1:N)'),cos(.47*(1:N)'));
verifyEqual(tc,O([],'size'),[N N]);
verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X),1e-10*norm(X));
[A,p]=cluster_floquet();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'scalable')));
verifyError(tc,@()cluster_floquet(5,.4,'periodic'), ...
    'cluster_floquet:Boundary');
verifyError(tc,@()cluster_floquet(4,ones(3,1)), ...
    'cluster_floquet:Angles');
verifyError(tc,@()cluster_floquet(4,.3,'invalid'), ...
    'cluster_floquet:Boundary');
verifyError(tc,@()cluster_floquet(4,.3,'open',.1,'bad'), ...
    'cluster_floquet:Mode');
verifyError(tc,@()cluster_floquet(11), ...
    'cluster_floquet:ExplicitOrder');
verifyError(tc,@()O(ones(N-1,1),'notransp'), ...
    'um_circuit_apply:Dimension');
verifyError(tc,@()O(ones(N,1),'bogus'), ...
    'um_circuit_apply:Flag');
end

function [Go,Ge]=symmetries(n)
odd=repmat({'I'},1,n);
even=odd;
odd(1:2:end)=repmat({'X'},1,numel(1:2:n));
even(2:2:end)=repmat({'X'},1,numel(2:2:n));
Go=pauli_word(odd);
Ge=pauli_word(even);
end

function U=pauli_word(names)
Z=[1 0;0 -1];X=[0 1;1 0];
U=1;
for j=1:numel(names)
    switch names{j}
        case 'I',G=eye(2);
        case 'X',G=X;
        case 'Z',G=Z;
        otherwise,error('test_cluster_floquet:Pauli','Unknown Pauli');
    end
    U=kron(U,G);
end
end
