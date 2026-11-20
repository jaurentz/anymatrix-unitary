function tests=test_circle_polynomial_unitary
tests=functiontests(localfunctions);
end

function test_polynomial_and_circle_spectrum(tc)
for n=[1 2 3 5 8 12]
    phases=2*pi*((0:n-1)'+.17)/n+.04*cos((1:n)');
    z=exp(1i*phases);
    c=poly(z);
    w=(1:n)';w=w/sum(w);
    for rep={'hessenberg','cmv'}
        name=rep{1};
        U=full(circle_polynomial_unitary(n,c,w,name));
        verifyEqual(tc,size(U),[n n]);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),8e-10*n);
        if strcmp(name,'hessenberg')
            verifyLessThan(tc,norm(tril(U,-2),'fro'),8e-11*n);
        end
        eigs=eig(U);
        verifyLessThan(tc,max(abs(polyval(c,eigs))),1e-8*n);
        [V,D]=eig(U);lambda=diag(D);
        for j=1:n
            [dist,k]=min(abs(lambda-z(j)));
            verifyLessThan(tc,dist,7e-9);
            verifyLessThan(tc,abs(abs(V(1,k))^2-w(j)),6e-8);
        end
        params=circle_polynomial_unitary(n,c,w,name,'parameters');
        verifyLessThan(tc,norm(params.polynomial-c),4e-10*n);
        verifyLessThan(tc,abs(params.beta- ...
            (-1)^(n+1)*conj(prod(z))),5e-10*n);
        F=circle_polynomial_unitary(n,c,w,name,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),5e-10*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),5e-10*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),5e-10*n);
    end
end
end

function test_strict_pencil_and_measure_separation(tc)
n=5;
z=exp(1i*([.1 1.2 2.4 3.7 5.1]'));
c=poly(z);w=(1:n)';w=w/sum(w);
S=circle_polynomial_unitary(n,c,w,'hessenberg','spectral');
verifyEqual(tc,S.finite.kind,'exact');
verifyEqual(tc,S.cyclic.kind,'exact');
verifyGreaterThan(tc,abs(S.cyclic.moments(1)-S.finite.moments(1)),.02);
P=circle_polynomial_unitary(n,c,w,'hessenberg','pencil');
verifyEqual(tc,size(P.A),[n n]);
verifyEqual(tc,size(P.B),[n n]);
eigp=eig(P.A,P.B);
for j=1:n
    verifyLessThan(tc,min(abs(eigp-z(j))),1e-8);
end
verifyLessThan(tc,norm((P.B\P.A)'*(P.B\P.A)-eye(n),'fro'),3e-10);
end

function test_self_inversive_not_sufficient_and_monomial_companion(tc)
% z^2-3z+1 is reciprocal but has non-unimodular roots.
verifyError(tc,@()circle_polynomial_unitary(2,[1 -3 1]), ...
    'circle_polynomial_unitary:CircleRoots');
verifyError(tc,@()circle_polynomial_unitary(2,[1 -2 1]), ...
    'circle_polynomial_unitary:Multiplicity');
verifyError(tc,@()circle_polynomial_unitary(2,[1 .2 2]), ...
    'circle_polynomial_unitary:SelfInversive');
z=exp(1i*[.15 1.2 2.8]);
c=poly(z);
C=[zeros(2,1),eye(2);-fliplr(c(2:end))];
verifyGreaterThan(tc,norm(C'*C-eye(3),'fro'),.1);
U=circle_polynomial_unitary(3,c);
verifyLessThan(tc,norm(full(U'*U-eye(3)),'fro'),4e-12);
end

function test_bad_input_and_defaults(tc)
[A,props]=circle_polynomial_unitary();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyEqual(tc,size(circle_polynomial_unitary(6)),[6 6]);
verifyError(tc,@()circle_polynomial_unitary(3,[1 2]), ...
    'circle_polynomial_unitary:Polynomial');
verifyError(tc,@()circle_polynomial_unitary(3,[1 0 0 0]), ...
    'circle_polynomial_unitary:SelfInversive');
verifyError(tc,@()circle_polynomial_unitary(3,[],[],'bad'), ...
    'circle_polynomial_unitary:Representation');
verifyError(tc,@()circle_polynomial_unitary(3,[],[],'cmv','bad'), ...
    'circle_polynomial_unitary:Mode');
end
