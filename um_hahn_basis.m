function [V,J,lambda]=um_hahn_basis(n,alpha,beta)
%UM_HAHN_BASIS Orthonormal Hahn functions via their self-adjoint Jacobi
%   difference operator. V has row index x and column degree k.
%   Q_k(x;alpha,beta,N)=3F2(-k,k+alpha+beta+1,-x;alpha+1,-N;1),
%   and the sign convention Q_k(0)=1 fixes V(1,k)>0.
%   The normalized dual-Hahn evaluation matrix is V' (dual roles).
%   https://doi.org/10.1007/978-3-642-05014-5
N=n-1;
if n==1
    V=1;J=0;lambda=0;
    return
end
x=(0:N)';
B=(x+alpha+1).*(x-N);
D=x.*(x-beta-N-1);
b=-sqrt(B(1:N).*D(2:end)); % reversible detailed-balance offdiagonal
J=diag(-B-D)+diag(b,1)+diag(b,-1);
[V,L]=eig(J,'vector');
[lambda,ord]=sort(real(L));
V=real(V(:,ord));
sgn=sign(V(1,:));
if any(sgn==0)
    error('um_hahn_basis:Precision', ...
        'Unable to fix the normalized Hahn polynomial signs.');
end
V=bsxfun(@times,V,sgn);
end
