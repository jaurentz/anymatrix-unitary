function tests=test_rational_arnoldi_circle
%TEST_RATIONAL_ARNOLDI_CIRCLE Finite orthogonal rational Krylov closure.
tests=functiontests(localfunctions);
end

function test_exact_full_unitarity_and_independent_spectral_data(tc)
for n=[1 2 4 7 12 19]
    z=exp(2i*pi*((0:n-1)'+.17)/n);
    w=(1:n)'.^.7;w=w/sum(w);
    families={inf(n-1,1),zeros(n-1,1), ...
        (2+.25i)*ones(n-1,1),(.3+.1i)*ones(n-1,1)};
    for f=1:numel(families)
        poles=families{f};
        T=rational_arnoldi_circle(n,poles,z,w);
        P=rational_arnoldi_circle(n,poles,z,w,'parameters');
        S=rational_arnoldi_circle(n,poles,z,w,'spectral');
        verifyEqual(tc,size(T),[n n]);
        verifyLessThan(tc,norm(T'*T-eye(n),'fro'),3e-10*n);
        verifyLessThan(tc,norm(P.Q'*P.Q-eye(n),'fro'),3e-10*n);
        verifyLessThan(tc,norm(P.unitary_closure-T,'fro'),4e-13*n);
        verifyLessThan(tc,norm(P.Q'*diag(z)*P.Q-T,'fro'),5e-12*n);
        verifyLessThan(tc,abs(sum(P.weights)-1),2e-13);
        verifyLessThan(tc,norm(P.Q(:,1)-sqrt(w)),2e-14*n);
        verifyLessThan(tc,norm(P.weighted_rational_basis'* ...
            diag(w)*P.weighted_rational_basis-eye(n),'fro'), ...
            3e-10*n);
        verifyEqual(tc,S.family,'rational_arnoldi_circle');
        verifyEqual(tc,S.finite.kind,'exact');
        verifyLessThan(tc,norm(S.finite.eigenvalues()-z),2e-13*n);
        verifyLessThan(tc,norm(S.cyclic.weights-w),2e-13*n);
        for k=-3:3
            verifyLessThan(tc,abs(S.finite.moments(k)- ...
                trace(T^k)/n),8e-9*n);
            verifyLessThan(tc,abs(S.cyclic.moments(k)- ...
                sum(w.*z.^k)),1e-12*n);
        end
        lam=eig(T);
        for j=1:n
            verifyLessThan(tc,min(abs(lam-z(j))),3e-8*n);
        end
    end
end
end

function test_rectangular_hessenberg_pencil(tc)
n=14;z=exp(2i*pi*((0:n-1)'+.11)/n);
w=(1:n)'.^1.4;w=w/sum(w);
cases={inf(n-1,1),zeros(n-1,1), ...
    (2+.3i)*ones(n-1,1), ...
    [Inf;2+1i;.4+.2i;3-.5i;-.2;.7i;-3; ...
     2.1+.2i;.3-.5i;Inf;3.1;-.1i;1.9]};
for j=1:numel(cases)
    poles=cases{j};
    B=rational_arnoldi_circle(n,poles,z,w,'basis');
    P=rational_arnoldi_circle(n,poles,z,w,'pencil');
    verifyEqual(tc,size(P.H),[n n-1]);
    verifyEqual(tc,size(P.K),[n n-1]);
    verifyFalse(tc,B.pencil_is_square);
    verifyLessThan(tc,norm(P.H-B.H,'fro'),1e-14*n);
    verifyLessThan(tc,norm(P.K-B.K,'fro'),1e-14*n);
    verifyLessThan(tc,norm(tril(P.H,-2),'fro'),4e-12*n);
    verifyLessThan(tc,norm(tril(P.K,-2),'fro'),4e-12*n);
    verifyLessThan(tc,norm(B.unitary_closure*P.K-P.H,'fro'), ...
        8e-10*n);
    verifyLessThan(tc,norm(diag(z)*B.Q*P.K-B.Q*P.H,'fro'), ...
        8e-10*n);
    verifyLessThan(tc,P.recurrence_residual,8e-10*n);
    if all(isinf(poles))
        verifyLessThan(tc,norm(tril(B.unitary_closure,-2),'fro'), ...
            6e-10*n);
        verifyLessThan(tc,norm(P.K-eye(n,n-1),'fro'),3e-12*n);
    end
end
end

function test_poles_nearly_coincident_and_reciprocal(tc)
n=10;z=exp(2i*pi*((0:n-1)'+.27)/n);
w=(1:n)';w=w/sum(w);
p=2+.45i;q=1/conj(p);
poles=[p;q;p+1e-8;q+2e-8;Inf; ...
    2.0000001+.45i;.4+.1i;3-.25i;-.15];
B=rational_arnoldi_circle(n,poles,z,w,'basis');
verifyLessThan(tc,norm(B.Q'*B.Q-eye(n),'fro'),1e-10*n);
verifyLessThan(tc,norm(B.unitary_closure'* ...
    B.unitary_closure-eye(n),'fro'),1e-10*n);
verifyLessThan(tc,B.recurrence_residual,2e-10*n);
verifyEqual(tc,B.poles,poles);
end

function test_polynomial_and_rational_models_distinct(tc)
n=9;z=exp(2i*pi*((0:n-1)'+.15)/n);
w=(1:n)';w=w/sum(w);
A=rational_arnoldi_circle(n,inf(n-1,1),z,w);
B=rational_arnoldi_circle(n,2*ones(n-1,1),z,w);
verifyGreaterThan(tc,norm(A-B,'fro'),.2);
% Same exact finite spectrum and initial cyclic measure for both.
S=rational_arnoldi_circle(n,2*ones(n-1,1),z,w,'spectral');
verifyLessThan(tc,abs(S.cyclic.moments(2)- ...
    sum(w.*z.^2)),2e-12);
end

function test_failures(tc)
n=5;z=exp(2i*pi*((0:n-1)'+.17)/n);
verifyError(tc,@()rational_arnoldi_circle(n,ones(n-1,1)), ...
    'rational_arnoldi_circle:Poles');
verifyError(tc,@()rational_arnoldi_circle(n,inf(n-1,1), ...
    [z(1);z(1);z(3:end)]), ...
    'rational_arnoldi_circle:Distinct');
verifyError(tc,@()rational_arnoldi_circle(n,inf(n-1,1),z, ...
    [1;0;2;3;4]),'rational_arnoldi_circle:Weights');
verifyError(tc,@()rational_arnoldi_circle(n,inf(n-1,1),z,[], ...
    'operator'),'rational_arnoldi_circle:Mode');
% A negative Laurent monomial followed by multiplication by z
% returns to the initial constant vector and causes exact breakdown.
verifyError(tc,@()rational_arnoldi_circle(n, ...
    [0;Inf;Inf;Inf],z,ones(n,1)), ...
    'rational_arnoldi_circle:Breakdown');
end
