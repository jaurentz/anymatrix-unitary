function [U,properties]=associated_legendre_transform(n,m)
%ASSOCIATED_LEGENDRE_TRANSFORM Fixed-order orthogonal Legendre block.
%   U=ASSOCIATED_LEGENDRE_TRANSFORM(N,M), integer M>=0, is the N-by-N
%   Gauss-Jacobi evaluation transform for Jacobi weight (1-x^2)^M.
%   Equivalently, the weighted associated-Legendre block consists
%   of degrees ell=M,...,M+N-1 at the corresponding N Gauss-Jacobi
%   nodes; the Condon-Shortley global sign (-1)^M is omitted.
%   This square fixed-order block is exactly orthogonal. It is NOT
%   the rectangular all-orders spherical-harmonic sampling transform.
%   The Jacobi matrix has zero diagonal and offdiagonal
%       sqrt(k*(k+2*M)/((2*k+2*M)^2-1)), k=1,...,N-1.
%   O(N^3) setup via eigenvectors; fast SHT kernels are not included.
%   References: Golub-Welsch (1969); Slevinsky (2017); Olver et al. (2020).
%   [93] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [95] https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D
%   [96] https://arxiv.org/abs/1705.05448
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(m),m=0;end
validateattributes(m,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
k=(1:n-1)';
b=sqrt((k.*(k+2*m))./((2*k+2*m).^2-1));
U=um_gauss_basis(zeros(n,1),b);
end
