function tests=test_fermionic_exterior_power
tests=functiontests(localfunctions);
end
function test_endpoint_and_one_particle(tc)
for d=[1 2 3 5]
    V=discrete_fourier(d);
    verifyEqual(tc,fermionic_exterior_power(d,0,V),1,'AbsTol',1e-13);
    verifyLessThan(tc,norm(fermionic_exterior_power(d,1,V)-V,'fro'),2e-13*d);
    verifyLessThan(tc,abs(fermionic_exterior_power(d,d,V)-det(V)), ...
        2e-13*d);
    for p=0:d
        U=fermionic_exterior_power(d,p,V);
        n=nchoosek(d,p);
        verifyEqual(tc,size(U),[n n]);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),5e-11*n);
        S=fermionic_exterior_power(d,p,V,'spectral');
        verifyEqual(tc,S.finite.kind,'exact');
        z=S.finite.eigenvalues();computed=eig(U);
        verifyEqual(tc,numel(z),n);
        for j=1:n
            verifyLessThan(tc,min(abs(computed(j)-z)),5e-10);
        end
        for k=-4:5
            verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/n), ...
                3e-10);
        end
    end
end
end
function test_group_representation_and_cofactor_signs(tc)
d=4;V=discrete_fourier(d);
W=diag(exp(1i*[.1 -.4 .25 .7]));
for p=0:d
    A=fermionic_exterior_power(d,p,V);
    B=fermionic_exterior_power(d,p,W);
    C=fermionic_exterior_power(d,p,V*W);
    verifyLessThan(tc,norm(A*B-C,'fro'),3e-12*d);
end
p=d-1;A=fermionic_exterior_power(d,p,V);
sets=nchoosek(1:d,p);
missing=zeros(d,1);
for k=1:d,missing(k)=setdiff(1:d,sets(k,:));end
for a=1:d
    for b=1:d
        i=missing(a);j=missing(b);
        expected=(-1)^(i+j)*det(V)*conj(V(i,j));
        verifyLessThan(tc,abs(A(a,b)-expected),3e-12);
    end
end
end
function test_distinct_occupancy_and_large_moments(tc)
d=4;z=exp(1i*[.15;.35;.8;1.2]);V=diag(z);
S=fermionic_exterior_power(d,2,V,'spectral');
idx=nchoosek(1:d,2);
expected=prod(z(idx),2);
verifyLessThan(tc,max(abs(sort(angle(S.finite.eigenvalues()))- ...
    sort(angle(expected)))),2e-12);
Sbig=fermionic_exterior_power(25,12,eye(25),'spectral');
verifyEqual(tc,Sbig.dimension,nchoosek(25,12));
verifyLessThan(tc,abs(Sbig.finite.moments(0)-1),1e-12);
verifyLessThan(tc,abs(Sbig.finite.moments(5)-1),1e-12);
end
function test_validation(tc)
verifyError(tc,@()fermionic_exterior_power(3,4), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()fermionic_exterior_power(3,2,ones(3,3)), ...
    'fermionic_exterior_power:Unitary');
verifyError(tc,@()fermionic_exterior_power(3,2,eye(2)), ...
    'fermionic_exterior_power:Base');
verifyError(tc,@()fermionic_exterior_power(3,2,eye(3),'operator'), ...
    'fermionic_exterior_power:Mode');
verifyError(tc,@()fermionic_exterior_power(12,6,eye(12)), ...
    'fermionic_exterior_power:Order');
end
