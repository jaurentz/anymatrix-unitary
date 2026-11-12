function [nodes,weights] = szego_quadrature(U)
%SZEGO_QUADRATURE Nodes and e1-weights of a finite unitary matrix.
%   [Z,W] = SZEGO_QUADRATURE(U) gives the eigenvalues Z and nonnegative
%   weights W=|v_j(1)|^2 from normalized eigenvectors. These weights
%   describe the distinguished vector measure, NOT the equal-weight
%   trace eigenvalue measure. Numerical eigendecomposition is used.
[nodes,weights]=js_szego_quadrature(U);
end
