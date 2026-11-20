function [U,properties]=gauss_hermite_transform(n)
%GAUSS_HERMITE_TRANSFORM Real orthogonal Gauss-Hermite evaluation matrix.
%   U=GAUSS_HERMITE_TRANSFORM(N) maps coefficients of the first N
%   orthonormal physicists' Hermite polynomials (weight exp(-x^2))
%   into quadrature-weighted function values at N Gauss-Hermite nodes.
%   Precisely U(j,k)=sqrt(w_j/sqrt(pi))*p_{k-1}(x_j) for polynomials
%   orthonormal under the normalized Gaussian weight exp(-x^2)/sqrt(pi).
%   U is exactly orthogonal in arithmetic over the reals.
%   The Jacobi matrix has diagonal zero and offdiagonal sqrt(k/2).
%   U's eigenphases are NOT the Gauss-Hermite nodes.
%   Uses a dense eigenvector construction; O(N^3) setup, O(N^2) multiply.
%   References: Golub and Welsch (1969); Webb (2026), arXiv:2604.02041.
%   [93] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [98] https://arxiv.org/abs/2604.02041
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
U=um_gauss_basis(zeros(n,1),sqrt((1:n-1)'/2));
end
