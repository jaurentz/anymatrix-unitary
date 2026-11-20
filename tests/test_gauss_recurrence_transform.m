function tests=test_gauss_recurrence_transform
tests=functiontests(localfunctions);
end
function test_general_positive_recursions(tc)
for n=[1 2 3 8 22]
    b=.25*sin((1:n)')+.1*(1:n)'/n;
    a=.3+.08*(1:n-1)';
    Q=gauss_recurrence_transform(n,b,a);
    verifyEqual(tc,size(Q),[n n]);
    verifyTrue(tc,isreal(Q));
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),2e-11*n);
    verifyGreaterThanOrEqual(tc,min(Q(:,1)),-3e-14);
    J=diag(b);
    if n>1,J=J+diag(a,1)+diag(a,-1);end
    lambda=sort(eig(J));
    verifyLessThan(tc,norm(Q*J*Q'-diag(lambda),'fro'),3e-11*n);
    verifyLessThan(tc,norm(Q-gauss_recurrence_transform(n,b,a),'fro'),2e-12*n);
end
end
function test_classical_existing_specializations(tc)
for n=[1 2 8 19]
    h=gauss_recurrence_transform(n);
    verifyLessThan(tc,norm(h-gauss_hermite_transform(n),'fro'),2e-12*n);
    a=sqrt((1:n-1)'/2);
    verifyLessThan(tc,norm(h-gauss_recurrence_transform(n,zeros(n,1),a),'fro'), ...
        2e-12*n);
    ell=1.4;k=(1:n-1)';
    b=(2*(0:n-1)'+ell+1);
    a=sqrt(k.*(k+ell));
    L=gauss_recurrence_transform(n,b,a);
    verifyLessThan(tc,norm(L-gauss_laguerre_transform(n,ell),'fro'),2e-11*n);
end
end
function test_unknown_unitary_spectral_law(tc)
S=gauss_recurrence_transform(10,[1:10]',ones(9,1),'spectral');
verifyEqual(tc,S.dimension,10);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyTrue(tc,contains(S.notes,'quadrature nodes'));
end
function test_input_validation(tc)
verifyError(tc,@()gauss_recurrence_transform(3,[1 2],ones(2,1)), ...
    'gauss_recurrence_transform:Diagonal');
verifyError(tc,@()gauss_recurrence_transform(3,[1 Inf 2],ones(2,1)), ...
    'gauss_recurrence_transform:Diagonal');
verifyError(tc,@()gauss_recurrence_transform(3,[1 2 3],[1 0]), ...
    'gauss_recurrence_transform:Offdiagonal');
verifyError(tc,@()gauss_recurrence_transform(3,[1 2 3],[-1 1]), ...
    'gauss_recurrence_transform:Offdiagonal');
verifyError(tc,@()gauss_recurrence_transform(3,[1 2 3],ones(3,1)), ...
    'gauss_recurrence_transform:Offdiagonal');
verifyError(tc,@()gauss_recurrence_transform(4,[],[],'operator'), ...
    'gauss_recurrence_transform:Mode');
end
