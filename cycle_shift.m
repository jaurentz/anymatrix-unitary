function [U, properties] = cycle_shift(n)
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
properties = {'unitary', 'orthogonal', 'real', 'permutation', ...
              'circulant', 'sparse', 'square', 'scalable', 'eigenvalues'};
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
j = 1:n;
i = mod(j,n) + 1;
U = sparse(i,j,1,n,n);
end
