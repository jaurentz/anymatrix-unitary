function tests=test_haar_tensor_certification
%TEST_HAAR_TENSOR_CERTIFICATION Tests for 2026 unitary family expansion.
tests=functiontests(localfunctions);
end

function test_interfaces_and_unitarity(testCase)
f={@haar_special_unitary,@haar_unimodular_coe,@haar_unimodular_cse, ...
   @entangled_reflection,@haar_tensor_power, ...
   @haar_block_tensor_power,@tensor_permutation};
for j=1:numel(f)
    [A,p]=f{j}();
    verifyEmpty(testCase,A);
    verifyTrue(testCase,any(strcmp(p,'unitary')));
    verifyTrue(testCase,any(strcmp(p,'square')));
end
for n=1:5
    checkunitary(testCase,haar_special_unitary(n,13));
    checkunitary(testCase,haar_unimodular_coe(n,13));
    checkunitary(testCase,haar_unimodular_cse(n,13));
    verifyLessThan(testCase,abs(det(haar_special_unitary(n,13))-1),1e-11);
    verifyLessThan(testCase,abs(det(haar_unimodular_coe(n,13))-1),1e-11);
    verifyLessThan(testCase,abs(det(haar_unimodular_cse(n,13))-1),1e-10);
end
end

function test_unimodular_symmetry(testCase)
for n=1:5
    A=haar_unimodular_coe(n,7);
    verifyLessThan(testCase,norm(A-A.','fro'),1e-11);
    B=haar_unimodular_cse(n,7);
    I=eye(n);
    J=[zeros(n),I;-I,zeros(n)];
    verifyLessThan(testCase,norm(B-J*B.'*J','fro'),2e-11);
    if n==1,verifyLessThan(testCase,norm(B-eye(2),'fro'),1e-11);end
end
end

function test_haar_seed_restoration(testCase)
state=rng;
A=haar_special_unitary(4,72);
B=haar_unimodular_coe(4,72);
C=haar_unimodular_cse(4,72);
D=haar_tensor_power(2,3,72);
E=haar_block_tensor_power(2,1,2,72);
verifyEqual(testCase,state,rng);
verifyEqual(testCase,A,haar_special_unitary(4,72),'AbsTol',0);
verifyEqual(testCase,B,haar_unimodular_coe(4,72),'AbsTol',0);
verifyEqual(testCase,C,haar_unimodular_cse(4,72),'AbsTol',0);
verifyEqual(testCase,D,haar_tensor_power(2,3,72),'AbsTol',0);
verifyEqual(testCase,E,haar_block_tensor_power(2,1,2,72),'AbsTol',0);
end

function test_entangled_reflection(testCase)
for d=2:6
    A=entangled_reflection(d);
    verifyTrue(testCase,issparse(A));
    checkunitary(testCase,A);
    verifyLessThan(testCase,norm(A-A','fro'),1e-12);
    ev=sort(real(eig(full(A))));
    verifyLessThan(testCase,abs(ev(1)+1),1e-12);
    verifyLessThan(testCase,norm(ev(2:end)-1),1e-12);
    O=entangled_reflection(d,pi,'operator');
    verifyEqual(testCase,O([], 'size'),[d*d d*d]);
    X=randn(d*d,3)+1i*randn(d*d,3);
    verifyLessThan(testCase,norm(O(X,'notransp')-A*X,'fro'),1e-12);
    verifyLessThan(testCase,norm(O(X,'ctransp')-A'*X,'fro'),1e-12);
end
A=entangled_reflection(4,.47);
checkunitary(testCase,A);
v=zeros(16,1);v(1:5:16)=1/2;
verifyLessThan(testCase,norm(A*v-exp(1i*.47)*v),1e-12);
end

function test_tensor_spectral_moments(testCase)
V=haar_unitary(3,24);
U=haar_tensor_power(3,3,24);
checkunitary(testCase,U);
for k=1:4
    verifyLessThan(testCase,abs(trace(U^k)/27-(trace(V^k)/3)^3),2e-12);
end
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(24,'twister');
V1=haar_unitary(2);V2=haar_unitary(2);
W=haar_block_tensor_power(2,1,2,24);
verifyLessThan(testCase,norm(W-kron(V1,kron(V2,V2)),'fro'),2e-12);
end

function test_tensor_permutation_and_trace(testCase)
for d=[2 3]
    S=tensor_permutation(d,[2 1]);
    checkunitary(testCase,S);
    verifyTrue(testCase,isreal(S));
    verifyTrue(testCase,issparse(S));
    verifyLessThan(testCase,norm(S*S-eye(d*d),'fro'),1e-12);
    verifyEqual(testCase,trace(S),d);
    verifyEqual(testCase,trace(S^2),d^2);
end
P=tensor_permutation(2,[2 3 1]);
verifyEqual(testCase,trace(P),2);
verifyEqual(testCase,trace(P^2),2);
verifyEqual(testCase,trace(P^3),8);
checkunitary(testCase,P);
verifyError(testCase,@()tensor_permutation(2,[1 1]), ...
    'tensor_permutation:Permutation');
end

function checkunitary(testCase,A)
n=size(A,1);
verifyEqual(testCase,size(A),[n n]);
verifyLessThan(testCase,norm(full(A'*A-eye(n)),'fro'),3e-11*n);
end
