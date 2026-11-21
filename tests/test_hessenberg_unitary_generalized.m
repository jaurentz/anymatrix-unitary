function tests=test_hessenberg_unitary_generalized
%TEST_HESSENBERG_UNITARY_GENERALIZED Nonunitary RQR benchmark separation.
tests=functiontests(localfunctions);
end

function test_nonunitary_quotient_and_generalized_spectrum(tc)
for n=[3 7 11 19]
    j=(1:n-1)';
    poles=(.45+.12*j).*exp(.17i*j);
    strict=unitary_hessenberg_pencil(n,poles,4,'complex','pencil');
    P=um_hessenberg_unitary_generalized(n,poles,4,.37);
    verifyEqual(tc,P.classification, ...
        'generalized_hessenberg_unitary_hessenberg');
    verifyFalse(tc,P.strict_unitary_pencil);
    verifyFalse(tc,P.unitary_factors);
    verifyLessThan(tc,norm(P.B'*P.B-eye(n),'fro'), ...
        2e-11*n);
    verifyLessThan(tc,norm(tril(P.A,-2),'fro'),2e-12*n);
    verifyLessThan(tc,norm(tril(P.B,-2),'fro'),2e-12*n);
    verifyLessThan(tc,norm(P.A-strict.A- ...
        .37*diag((1:n)'/(n+1)),'fro'),2e-12*n);
    verifyLessThan(tc,norm(P.B-strict.B,'fro'),2e-12*n);
    verifyLessThan(tc,norm(P.U-P.B'*P.A,'fro'),2e-12*n);
    verifyGreaterThan(tc,norm(P.A'*P.A-eye(n),'fro'),1e-3);
    verifyGreaterThan(tc,norm(P.U'*P.U-eye(n),'fro'),1e-3);
    verifyLessThan(tc,norm(P.poles-poles),2e-9*n);
    z1=eig(P.A,P.B);z2=eig(P.U);
    for k=1:n
        verifyLessThan(tc,min(abs(z1-z2(k))),2e-9*n);
    end
end
end

function test_real_nonunitary_pencil_with_projective_poles(tc)
poles=[-2;0;Inf;1.6;.2;-Inf];
n=numel(poles)+1;
P=um_hessenberg_unitary_generalized(n,poles,20,.5,'real');
verifyTrue(tc,isreal(P.A));
verifyTrue(tc,isreal(P.B));
verifyTrue(tc,isreal(P.U));
verifyEqual(tc,P.poles(2),0);
verifyTrue(tc,isinf(P.poles(3)));
verifyTrue(tc,isinf(P.poles(6)));
verifyLessThan(tc,norm(P.B'*P.B-eye(n),'fro'),2e-11*n);
verifyGreaterThan(tc,norm(P.U'*P.U-eye(n),'fro'),1e-3);
lambda=eig(P.U);
for j=1:n
    verifyLessThan(tc,min(abs(lambda-conj(lambda(j)))),5e-10*n);
end
verifyError(tc,@()um_hessenberg_unitary_generalized(n,poles,20,0), ...
    'MATLAB:expectedPositive');
end
