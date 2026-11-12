function [nodes,weights] = js_szego_quadrature(U)
%JS_SZEGO_QUADRATURE Quadrature nodes and e1 spectral weights of a unitary.
%   [NODES,WEIGHTS]=JS_SZEGO_QUADRATURE(U) computes eigenvalues and the
%   squared absolute value of the first entry of normalized eigenvectors.
%   For an OPUC finite closure these integrate Laurent polynomials of
%   degree -(n-1),...,(n-1) exactly (up to floating-point error).
validateattributes(U,{'numeric'},{'2d','square'});
[V,D]=eig(full(U));
nodes=diag(D);
weights=abs(V(1,:).').^2;
[~,idx]=sort(angle(nodes));
nodes=nodes(idx);
weights=weights(idx);
end
