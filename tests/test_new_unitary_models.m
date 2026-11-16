function tests=test_new_unitary_models
%TEST_NEW_UNITARY_MODELS Independent scattering and magnetic-walk checks.
tests=functiontests(localfunctions);
end

function test_scattering_default_and_pencil(testCase)
for n=[2 4 8]
    for q=[1 2 3]
        Z=scattering_zipper(n,q);
        P=scattering_zipper(n,q,[],[],[],'pencil');
        verifyEqual(testCase,size(Z),[n*q n*q]);
        verifyTrue(testCase,issparse(Z));
        verifyLessThan(testCase,norm(full(Z'*Z-speye(n*q)),'fro'),5e-11);
        verifyLessThan(testCase,norm(full(P.B\P.A-Z),'fro'),5e-11);
    end
end
end

function test_scattering_direct_layers(testCase)
q=2;n=4;
a=[.25 .1i;-.05 .13];
r=sqrtm(eye(q)-a*a');rt=sqrtm(eye(q)-a'*a);
F=[1 1;1 -1]/sqrt(2);
S=[a,r*F;F'*rt,-F'*a'*F];
blocks=repmat(S,1,1,n-1);
L=diag(exp(1i*[.1 .3]));
R=[0 1;-1 0];
P=scattering_zipper(n,q,blocks,L,R,'pencil');
V=blkdiag(S,S);
W=blkdiag(L,S,R);
verifyLessThan(testCase,norm(full(P.U-V*W),'fro'),2e-12);
verifyLessThan(testCase,norm(full(P.A-W),'fro'),2e-12);
verifyLessThan(testCase,norm(full(P.B-V'),'fro'),2e-12);
end

function test_scattering_exact_cycle_limit(testCase)
for n=[2 4 8]
    for q=[1 2 3]
        S=[zeros(q),eye(q);eye(q),zeros(q)];
        Z=scattering_zipper(n,q,repmat(S,1,1,n-1));
        verifyLessThan(testCase,norm(full(Z^n-speye(n*q)),'fro'),1e-12);
        for k=1:n-1
            verifyLessThan(testCase,abs(trace(Z^k)),1e-12);
        end
    end
end
end

function test_uamo_unitarity_and_pencil(testCase)
for n=[1 2 3 5 9]
    for pair=[0 0;1 0;1 1;.2 .7;.7 .2;.67 .67]'
        W=unitary_almost_mathieu(n,pair(1),pair(2),2/5,.23);
        P=unitary_almost_mathieu(n,pair(1),pair(2),2/5,.23,'pencil');
        verifyTrue(testCase,issparse(W));
        verifyEqual(testCase,size(W),[2*n 2*n]);
        verifyLessThan(testCase,norm(full(W'*W-speye(2*n)),'fro'), ...
            3e-12*n);
        verifyLessThan(testCase,norm(full(P.B\P.A-W),'fro'), ...
            3e-12*n);
        verifyLessThan(testCase,norm(full(P.A'*P.A-speye(2*n)),'fro'), ...
            3e-12*n);
    end
end
end

function test_uamo_direct_action(testCase)
n=7;l1=.63;l2=.91;freq=(sqrt(5)-1)/2;phase=.18;
x=(1:2*n)'+1i*(2*n:-1:1)';
actual=unitary_almost_mathieu(n,l1,l2,freq,phase)*x;
states=reshape(x,2,n);
t=2*pi*((0:n-1)*freq+phase);
r1=sqrt(1-l1^2);r2=sqrt(1-l2^2);
qp=(l2*cos(t)+1i*r2).*states(1,:)-l2*sin(t).*states(2,:);
qm=l2*sin(t).*states(1,:)+(l2*cos(t)-1i*r2).*states(2,:);
expected=zeros(2,n);
for j=1:n
    prev=mod(j-2,n)+1;
    next=mod(j,n)+1;
    expected(1,j)=l1*qp(prev)-r1*qm(j);
    expected(2,j)=l1*qm(next)+r1*qp(j);
end
verifyLessThan(testCase,norm(actual-expected(:)),1e-12);
end

function test_uamo_free_shift_spectrum(testCase)
for n=[1 2 5 9]
    W=unitary_almost_mathieu(n,1,0,.367,.26);
    roots=exp(2i*pi*(0:n-1)'/n);
    exact=[1i*roots;-1i*conj(roots)];
    observed=eig(full(W));
    verifyLessThan(testCase,max(arrayfun(@(z)min(abs(z-exact)),observed)), ...
        3e-12);
end
end

function test_descriptor_contracts(testCase)
S=scattering_zipper(12,2,[],[],[],'spectral');
verifyEqual(testCase,S.family,'scattering_zipper');
verifyEqual(testCase,S.dimension,24);
verifyEqual(testCase,S.finite.kind,'unknown');
T=unitary_almost_mathieu(31,.7,.7,[],[],'spectral');
verifyEqual(testCase,T.family,'unitary_almost_mathieu');
verifyEqual(testCase,T.dimension,62);
verifyEqual(testCase,T.limit.kind,'unknown');
[e,p]=scattering_zipper();verifyEmpty(testCase,e);
verifyTrue(testCase,ismember('unitary',p));
[e,p]=unitary_almost_mathieu();verifyEmpty(testCase,e);
verifyTrue(testCase,ismember('unitary',p));
end

function test_parameter_validation(testCase)
verifyError(testCase,@()scattering_zipper(3,2), ...
            'scattering_zipper:EvenSites');
verifyError(testCase,@()scattering_zipper(4,2,ones(4,4,3)), ...
            'scattering_zipper:Blocks');
verifyError(testCase,@()scattering_zipper(4,2,[],2*eye(2)), ...
            'scattering_zipper:Boundary');
verifyError(testCase,@()unitary_almost_mathieu(4,.7,.5,.2,0,'operator'), ...
            'unitary_almost_mathieu:Mode');
end
