function [x,w]=um_legendre_gauss_rule(L)
%UM_LEGENDRE_GAUSS_RULE Ascending Gauss-Legendre abscissae and dx weights.
%   Exact-in-mathematics L-node quadrature for polynomials degree<=2L-1.
%   Sum(w)=2, x sorted ascending; no spectral phase interpretation.
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',64});
if L==1,x=0;w=2;return;end
k=(1:L-1)';
e=k./sqrt((2*k-1).*(2*k+1));
J=diag(e,1)+diag(e,-1);
[V,D]=eig(J);
[x,perm]=sort(real(diag(D)),'ascend');
w=2*real(V(1,perm)).^2;
x=x(:);w=w(:);
end
