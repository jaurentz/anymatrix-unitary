function [U, properties] = finite_cmv(n, alpha,output)
%FINITE_CMV Finite CMV multiplication matrix from Verblunsky coefficients.
%   U = FINITE_CMV(N, ALPHA) constructs an N-by-N CMV matrix, where
%   ALPHA is an N-vector of coefficients with |ALPHA(1:N-1)| < 1 and
%   |ALPHA(N)| = 1. The final coefficient is the finite boundary phase.
%   ALPHA defaults to [zeros(N-1,1); 1]. For N=1, U=conj(ALPHA(1)).
%   Uses the same shared core-construction machinery as JACOBI_SZEGO
%   and OPUC_MATRIX, with no duplicated core implementation.
%
%   Reference: Cantero, Moral, Velazquez, Linear Algebra Appl. 362 (2003).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [6] https://doi.org/10.1016/S0024-3795(02)00457-3
%   [7] https://doi.org/10.1016/j.cam.2006.10.033
properties = {'square','unitary','normal','sparse','banded','scalable', ...
              'parameter-dependent'};
if nargin == 0
    U = [];
    return;
end
validateattributes(n, {'numeric'}, {'scalar','integer','positive','finite'});
if nargin < 2 || isempty(alpha)
    alpha = [zeros(n-1,1); 1];
end
validateattributes(alpha, {'numeric'}, {'vector','finite'});
if numel(alpha) ~= n
    error('finite_cmv:DimensionMismatch', 'ALPHA must contain N coefficients.');
end
if nargin<3,output='matrix';end
U = js_core_product(alpha(1:end-1), alpha(end), 'cmv', [], output);
end
