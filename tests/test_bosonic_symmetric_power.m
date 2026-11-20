function tests=test_bosonic_symmetric_power
tests=functiontests(localfunctions);
end

function test_small_fock_unitarity_and_phases(tc)
for d=1:4
    V=discrete_fourier(d);
    for p=0:3
        U=bosonic_symmetric_power(d,p,V);
        N=nchoosek(d+p-1,p);
        verifyEqual(tc,size(U),[N N]);
        verifyLessThan(tc,norm(U'*U-eye(N),'fro'),8e-10*N);
        S=bosonic_symmetric_power(d,p,V,'spectral');
        verifyEqual(tc,S.dimension,N);
        verifyEqual(tc,S.finite.kind,'exact');
        z=S.finite.eigenvalues();
        verifyEqual(tc,numel(z),N);
        eigs=eig(U);
        for j=1:N
            verifyLessThan(tc,min(abs(eigs-z(j))),2e-8);
        end
        for k=-3:5
            verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/N), ...
                5e-9);
        end
    end
end
end

function test_creation_operator_factors_and_composition(tc)
d=2;p=2;V=discrete_fourier(d);
A=bosonic_symmetric_power(d,p,V);
% Occupations in separator order: [0,2], [1,1], [2,0].
verifyLessThan(tc,abs(A(2,2) - ...
    (V(1,1)*V(2,2)+V(1,2)*V(2,1))),2e-13);
verifyLessThan(tc,abs(A(3,2)-sqrt(2)*V(1,1)*V(1,2)),2e-13);
d=3;W=diag(exp(1i*[.2 .4 -.6]));V=discrete_fourier(d);
for p=0:3
    A=bosonic_symmetric_power(d,p,V);
    B=bosonic_symmetric_power(d,p,W);
    C=bosonic_symmetric_power(d,p,V*W);
    verifyLessThan(tc,norm(A*B-C,'fro'),1e-9);
end
end

function test_repeat_occupations_and_large_moments(tc)
z=exp(1i*[.12;.31;.87]);d=3;p=2;V=diag(z);
S=bosonic_symmetric_power(d,p,V,'spectral');
occ=um_boson_occupations(d,p);
expected=prod(bsxfun(@power,z(:).',occ),2);
actual=S.finite.eigenvalues();
verifyLessThan(tc,max(abs(sort(angle(actual))-sort(angle(expected)))),1e-12);
Sbig=bosonic_symmetric_power(20,12,eye(20),'spectral');
verifyEqual(tc,Sbig.dimension,nchoosek(31,12));
verifyLessThan(tc,abs(Sbig.finite.moments(0)-1),2e-12);
verifyLessThan(tc,abs(Sbig.finite.moments(5)-1),2e-12);
verifyError(tc,@()Sbig.finite.eigenvalues(), ...
    'um_boson_spectral:Materialization');
end

function test_validation_and_limits(tc)
verifyError(tc,@()bosonic_symmetric_power(3,2,ones(3)), ...
    'bosonic_symmetric_power:Unitary');
verifyError(tc,@()bosonic_symmetric_power(3,2,eye(2)), ...
    'bosonic_symmetric_power:Base');
verifyError(tc,@()bosonic_symmetric_power(12,6,eye(12)), ...
    'bosonic_symmetric_power:Order');
verifyError(tc,@()bosonic_symmetric_power(3,2,eye(3),'operator'), ...
    'bosonic_symmetric_power:Mode');
verifyEqual(tc,bosonic_symmetric_power(1,8,1i),(1i)^8,'AbsTol',1e-14);
end
