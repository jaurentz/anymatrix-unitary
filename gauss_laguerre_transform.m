function [U,properties]=gauss_laguerre_transform(n,alpha)
%GAUSS_LAGUERRE_TRANSFORM Orthogonal Gauss-Laguerre evaluation matrix.
%   U=GAUSS_LAGUERRE_TRANSFORM(N,ALPHA) uses generalized Laguerre
%   polynomials orthonormal with respect to x^alpha*exp(-x) on [0,inf).
%   ALPHA > -1; default ALPHA=0. Rows are Gauss nodes and columns
%   polynomial degrees 0:N-1. Quadrature weights are normalized to
%   sum one. The sign convention is positive leading coefficient,
%   rather than the usual (-1)^k leading coefficient of L_k^alpha.
%   The Jacobi diagonal is 2*k+alpha+1 (k=0:N-1); subdiagonal
%   sqrt(k*(k+alpha)) (k=1:N-1). Hence U'*U=I in exact arithmetic.
%   The eigenvalues of U are not the Laguerre quadrature nodes.
%   Explicit Golub-Welsch construction; no fast matvec is claimed.
%   Reference: Golub and Welsch, Math. Comp. 23 (1969), 221-230.
%   [93] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [95] https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=0;end
validateattributes(alpha,{'numeric'}, ...
    {'scalar','real','finite','>',-1});
k=(0:n-1)';
j=(1:n-1)';
U=um_gauss_basis(2*k+alpha+1,sqrt(j.*(j+alpha)));
end
