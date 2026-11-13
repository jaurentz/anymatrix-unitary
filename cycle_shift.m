function [U, properties] = cycle_shift(n, output)
%CYCLE_SHIFT Cyclic permutation matrix, a unitary test matrix.
%   U = CYCLE_SHIFT(N) returns the N-by-N cyclic right shift, defined by
%   U*e_j = e_{j+1} (indices are taken modulo N).
%
%   The eigenvalues are the Nth roots of unity. In particular, U^N = I,
%   and trace(U^k)/N is 1 if N divides k, and 0 otherwise.
%   U is stored in MATLAB's sparse format.
%
%   Reference:
%   P. J. Davis, Circulant Matrices, Wiley, 1979.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [1] https://openlibrary.org/books/OL4407364M/Circulant_matrices
properties = {'unitary', 'orthogonal', 'real', 'permutation', ...
              'circulant', 'sparse', 'square', 'scalable', 'eigenvalues'};
% Optional spectral-reference mode; preserve legacy output properties.
if nargin>=2 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('cycle_shift',n);
    return
end
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
j = 1:n;
i = mod(j,n) + 1;
if nargin<2,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_monomial_apply(X,flag,i,ones(n,1));
else
    U = sparse(i,j,1,n,n);
end
end
