function [U, properties] = diagonal_phase(n, theta, output)
%DIAGONAL_PHASE Sparse unitary diagonal matrix with prescribed eigenphases.
%   U = DIAGONAL_PHASE(N) has eigenvalues exp(2*pi*1i*k/N), k=0:N-1.
%   U = DIAGONAL_PHASE(N,THETA) uses the real N-vector of phases THETA.
%
%   The eigenvalues are exactly the diagonal entries, with their
%   multiplicities. Finite-dimensional spectra are discrete measures;
%   any continuous-looking density requires an explicit smoothing rule.
%
%   Uses: controlled clustered spectra, arc gaps, phase wrapping, and
%   exact Fourier moment benchmarks.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [1] https://openlibrary.org/books/OL4407364M/Circulant_matrices
properties = {'unitary', 'square', 'sparse', 'scalable', ...
              'complex', 'eigenvalues', 'parameter-dependent'};
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
if nargin < 2 || isempty(theta)
    theta = 2*pi*(0:n-1)'/n;
end
validateattributes(theta, {'numeric'}, {'vector', 'real', 'finite'});
if numel(theta) ~= n
    error('diagonal_phase:DimensionMismatch', ...
          'THETA must contain precisely N real entries.');
end
if nargin<3,output='matrix';end
phase=exp(1i*theta(:));
if um_is_operator(output)
    U=@(X,flag)um_monomial_apply(X,flag,1:n,phase);
else
    U=spdiags(phase,0,n,n);
end
end
