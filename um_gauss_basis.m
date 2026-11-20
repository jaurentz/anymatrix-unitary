function Q=um_gauss_basis(diagonal,offdiagonal)
%UM_GAUSS_BASIS Orthogonal Gauss quadrature evaluation via Golub-Welsch.
%   Q(j,k)=sqrt(w_j/mu0)*p_{k-1}(x_j), for orthonormal polynomials
%   relative to the normalized measure dmu/mu0. The Jacobi matrix
%   has positive offdiagonals (positive-leading-coefficient convention).
%   Q maps polynomial coefficients to quadrature-weighted values.
%   This is an explicit O(n^3) eigensolver construction, NOT a fast
%   transform implementation. It avoids underflow in quadrature weights.
%
%   G. H. Golub and J. H. Welsch, Math. Comp. 23 (1969), 221-230.
d=diagonal(:); e=offdiagonal(:);
n=numel(d);
assert(numel(e)==max(0,n-1),'Jacobi offdiagonal has incorrect length.');
if n==1
    Q=1;
    return
end
J=diag(d)+diag(e,1)+diag(e,-1);
[V,D]=eig(J);
[~,idx]=sort(diag(D),'ascend');
V=V(:,idx);
s=sign(V(1,:));
s(s==0)=1;
V=bsxfun(@times,V,s);
Q=V.';
end
