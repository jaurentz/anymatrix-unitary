function tests=test_lee_yang_ising
%TEST_LEE_YANG_ISING Independent Ising partition-function checks.
tests=functiontests(localfunctions);
end

function test_partition_zeros_and_unitarity(testCase)
% Independent 2^n spin enumeration, not the CMV/Szego recurrence.
allp=[.25 .61 .1 .43 .8 .28 .36]';
for n=1:7
    p=allp(1:n);
    U=lee_yang_ising(n,p,'custom');
    P=lee_yang_ising(n,p,'custom','pencil');
    verifyEqual(testCase,size(U),[2*n 2*n]);
    verifyTrue(testCase,issparse(U));
    verifyLessThan(testCase,norm(full(U'*U-speye(2*n)),'fro'), ...
        3e-12*n);
    verifyLessThan(testCase,norm(full(P.B\P.A-U),'fro'), ...
        3e-12*n);
    verifyLessThan(testCase,norm(full(P.A*P.A'-P.B*P.B'),'fro'), ...
        3e-12*n);
    c=zeros(n+1,1);
    for mask=0:2^n-1
        spins=2*double(bitget(uint32(mask),1:n))-1;
        neighbors=spins([2:n 1]);
        weight=exp(sum(p(:)'.*spins.*neighbors));
        idx=(n+sum(spins))/2+1;
        c(idx)=c(idx)+weight;
    end
    z=eig(full(U));
    for j=1:numel(z)
        % z^n * Z_n(z) = sum_{k=0}^n c_k*z^(2k).
        verifyLessThan(testCase,abs(polyval(flipud(c),z(j)^2)), ...
            1e-9*sum(c));
        verifyLessThan(testCase,min(abs(z+z(j))),3e-10);
    end
end
end

function test_fibonacci_and_other_parameterizations(testCase)
p=[2/3 1/100];
word=[1 2 1 1 2 1 2 1]; % S^4(a)=abaababa.
U=lee_yang_ising(8,p,'fibonacci');
V=lee_yang_ising(8,p(word),'custom');
verifyEqual(testCase,U,V,'AbsTol',0);
verifyLessThan(testCase,norm(full(U'*U-speye(16)),'fro'),1e-12);
a=lee_yang_ising(8,.4,'uniform');
b=lee_yang_ising(8,[.4,.4],'alternating');
c=lee_yang_ising(8,[.4,.4],'fibonacci');
verifyEqual(testCase,a,b,'AbsTol',0);
verifyEqual(testCase,a,c,'AbsTol',0);
[empty,properties]=lee_yang_ising();
verifyEmpty(testCase,empty);
verifyTrue(testCase,any(strcmp(properties,'unitary')));
verifyTrue(testCase,any(strcmp(properties,'sparse')));
end

function test_spectral_descriptor_and_input_validation(testCase)
S=lee_yang_ising(34,[],[],'spectral');
verifyEqual(testCase,S.family,'lee_yang_ising');
verifyEqual(testCase,S.dimension,68);
verifyEqual(testCase,S.finite.kind,'unknown');
verifyEqual(testCase,S.limit.kind,'unknown');
verifyError(testCase,@()lee_yang_ising(5,[-1,.3],'fibonacci'), ...
    'MATLAB:expectedPositive');
verifyError(testCase,@()lee_yang_ising(5,[.1,.2],'custom'), ...
    'lee_yang_ising:Couplings');
verifyError(testCase,@()lee_yang_ising(5,.4,'nonexistent'), ...
    'lee_yang_ising:Pattern');
verifyError(testCase,@()lee_yang_ising(5,.4,'uniform','operator'), ...
    'lee_yang_ising:Mode');
end
