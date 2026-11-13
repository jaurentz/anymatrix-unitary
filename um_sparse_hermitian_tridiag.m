function H=um_sparse_hermitian_tridiag(d,u)
%UM_SPARSE_HERMITIAN_TRIDIAG Sparse Hermitian tridiagonal from diagonals.
%   H=UM_SPARSE_HERMITIAN_TRIDIAG(D,U) has real diagonal D,
%   superdiagonal U, and subdiagonal conj(U). It stores O(n) entries.
%   EXPM(H) is generally dense; sparse H is only the generator.
d=d(:);
u=u(:);
n=numel(d);
if n<1 || numel(u)~=n-1
    error('um_sparse_hermitian_tridiag:Shape', ...
          'Expected n real diagonal and n-1 superdiagonal entries.');
end
if ~isreal(d) || any(~isfinite(d)) || any(~isfinite(u))
    error('um_sparse_hermitian_tridiag:Parameters', ...
          'The diagonal must be real and all entries finite.');
end
H=spdiags(d,0,n,n);
if n>1
    H=H+sparse(1:n-1,2:n,u,n,n)+ ...
        sparse(2:n,1:n-1,conj(u),n,n);
end
end
