function tests=test_clifford_circuit
%TEST_CLIFFORD_CIRCUIT Clifford normalizer, exact phases, implicit actions.
tests=functiontests(localfunctions);
end

function test_exact_presets_and_multiplicities(tc)
presets={'identity','hadamard','phase','phase_dagger', ...
    'pauli_x','pauli_y','pauli_z','cnot','cz','swap'};
for q=[1 2 3]
    for j=1:numel(presets)
        name=presets{j};
        if q==1 && j>=8,continue;end
        A=clifford_circuit(q,name);
        S=clifford_circuit(q,name,'spectral');
        N=2^q;
        verifyEqual(tc,size(A),[N N]);
        verifyLessThan(tc,norm(A'*A-eye(N),'fro'),3e-11);
        verifyEqual(tc,S.dimension,N);
        verifyEqual(tc,S.finite.kind,'exact');
        ref=S.finite.eigenvalues();
        verifyEqual(tc,numel(ref),N);
        lam=eig(A);
        z=exp(1i*S.finite.atoms.angles(:));
        w=S.finite.atoms.weights(:);
        for k=1:numel(z)
            verifyEqual(tc,sum(abs(lam-z(k))<1e-8),round(N*w(k)));
        end
        verifyEqual(tc,sum(w),1,'AbsTol',1e-14);
        for k=-3:5
            verifyLessThan(tc,abs(mean(lam.^k)-S.finite.moments(k)), ...
                5e-10);
        end
    end
end
end

function test_custom_word_embedded_gates_and_pauli_normalizer(tc)
q=3;
seq={'H',1;'S',2;'CNOT',[1 3];'CZ',[2 3]; ...
    'SDG',1;'X',2;'SWAP',[3 1];'Y',3;'Z',1};
A=clifford_circuit(q,seq);
verifyLessThan(tc,norm(A'*A-eye(8),'fro'),8e-12);
local={ [1 1;1 -1]/sqrt(2),diag([1 1i]), ...
    [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0], ...
    diag([1 1 1 -1]),diag([1 -1i]),[0 1;1 0], ...
    [1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1], ...
    [0 -1i;1i 0],diag([1 -1])};
B=eye(8);
for j=1:size(seq,1)
    sites=seq{j,2};G=local{j};
    if numel(sites)==2
        W=qm_two_site_gate(q,G,sites(1),sites(2));
    else
        W=1;
        for i=1:q
            if i==sites
                W=kron(W,G);
            else
                W=kron(W,eye(2));
            end
        end
    end
    B=W*B;
end
verifyLessThan(tc,norm(A-B,'fro'),8e-12);
paulis={eye(2),[0 1;1 0],[0 -1i;1i 0],diag([1 -1])};
words=cell(4^q,1);
for j=0:4^q-1
    W=1;
    for s=1:q
        digit=mod(floor(j/4^(q-s)),4)+1;
        W=kron(W,paulis{digit});
    end
    words{j+1}=W;
end
for site=1:q
    for p=2:4
        P=1;
        for j=1:q
            if j==site,P=kron(P,paulis{p});
            else,P=kron(P,eye(2));end
        end
        conjugated=A*P*A';
        coeff=zeros(numel(words),1);
        for j=1:numel(words)
            coeff(j)=trace(words{j}'*conjugated)/8;
        end
        verifyLessThan(tc,abs(max(abs(coeff))-1),5e-12);
        verifyLessThan(tc,max(abs(imag(coeff))),5e-12);
    end
end
S=clifford_circuit(q,seq,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
end

function test_operator_and_large_statevector(tc)
for q=[1 2 4 6]
    seq={'H',1;'S',q;'SDG',1;'Y',q};
    if q>=2,seq=[seq;{'CNOT',[1 q]};{'CZ',[q 1]}];end
    A=clifford_circuit(q,seq);
    O=clifford_circuit(q,seq,'operator');
    X=reshape(sin((1:3*2^q)*.7),2^q,3)+ ...
        1i*reshape(cos((1:3*2^q)*.4),2^q,3);
    verifyEqual(tc,O([],'size'),[2^q 2^q]);
    verifyLessThan(tc,norm(O(X,'notransp')-A*X,'fro'),1e-9);
    verifyLessThan(tc,norm(O(X,'ctransp')-A'*X,'fro'),1e-9);
    verifyLessThan(tc,norm(O(X,'transp')-A'*X,'fro'),1e-9);
    verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'),1e-9);
end
q=14; N=2^q;
seq={'H',1;'S',14;'CNOT',[1 14];'CZ',[6 10]};
O=clifford_circuit(q,seq,'operator');
X=complex(sin((1:N)'*.43),cos((1:N)'*.23));
verifyEqual(tc,O([],'size'),[N N]);
verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X),1e-10*norm(X));
end

function test_identities_and_input_validation(tc)
[A,props]=clifford_circuit();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyFalse(tc,any(strcmp(props,'orthogonal')));
verifyEqual(tc,clifford_circuit(2,'identity'),eye(4),'AbsTol',0);
verifyEqual(tc,clifford_circuit(1),[1 1;1 -1]/sqrt(2),'AbsTol',2e-14);
verifyEqual(tc,clifford_circuit(2,{}),eye(4),'AbsTol',0);
s=clifford_circuit(3,{'S',3},'spectral');
verifyEqual(tc,s.finite.kind,'exact');
verifyEqual(tc,sum(s.finite.atoms.weights),1,'AbsTol',0);
seq=struct('name',{'H','CNOT'},'sites',{2,[2 1]});
verifyLessThan(tc,norm(clifford_circuit(2,seq) ...
    -clifford_circuit(2,{'H',2;'CX',[2 1]}),'fro'),1e-13);
verifyError(tc,@()clifford_circuit(2,'bad'),'clifford_circuit:Circuit');
verifyError(tc,@()clifford_circuit(1,'cnot'),'clifford_circuit:Sites');
verifyError(tc,@()clifford_circuit(2,{'CNOT',[1 1]}), ...
    'clifford_circuit:Sites');
verifyError(tc,@()clifford_circuit(2,{'H',[1 2]}), ...
    'clifford_circuit:Sites');
verifyError(tc,@()clifford_circuit(2,{'T',1}), ...
    'clifford_circuit:Gate');
verifyError(tc,@()clifford_circuit(2,{'H',3}), ...
    'clifford_circuit:Sites');
verifyError(tc,@()clifford_circuit(11,'identity'), ...
    'clifford_circuit:ExplicitOrder');
verifyError(tc,@()clifford_circuit(2,'identity','bad'), ...
    'clifford_circuit:Mode');
O=clifford_circuit(2,'bell','operator');
verifyError(tc,@()O(ones(3,1),'notransp'),'um_circuit_apply:Dimension');
verifyError(tc,@()O(ones(4,1),'bad'),'um_circuit_apply:Flag');
end
