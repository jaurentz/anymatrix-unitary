function [Q,info]=um_gauss_basis(diagonal,offdiagonal)
%UM_GAUSS_BASIS Orthogonal Gauss transform via symmetric Jacobi eigenvectors.
%   Q(j,k)=sqrt(w_j/mu0)*p_(k-1)(x_j). The Jacobi matrix has
%   positive offdiagonals and the polynomials have positive leading
%   coefficient. Q maps polynomial coefficients to weighted samples.
%
%   For large Gaussian orders, the first eigenvector component can
%   UNDERFLOW to zero (especially near extreme Gauss-Hermite or
%   Gauss-Laguerre nodes). Orienting each eigenvector by sign(V(1,j))
%   then becomes unstable. The Sturm/interlacing sign relation
%       sign(V(n,j))=(-1)^(n-j)*sign(V(1,j))
%   for eigenvalues sorted ASCENDING (j=1,...,n) and strictly
%   positive Jacobi offdiagonals allows using the MORE MAGNITUDINALLY
%   RELIABLE endpoint to fix the polynomial sign instead.
%   Thus Q(:,1) is theoretically positive even when the tiny leading
%   quadrature weights cannot be represented in finite precision.
%
%   [Q,INFO]=UM_GAUSS_BASIS(B,A) also returns sorted quadrature nodes,
%   representable normalized weights Q(:,1).^2, independently
%   evaluated LOG normalized weights by a scaled three-term recurrence
%   (including true values below realmin), and Gram/diagonalization
%   residuals. Eigenphases of Q are NOT the real quadrature nodes.
%
%   This implementation still computes a complete dense eigensystem:
%   O(n^3) setup, O(n^2) storage, O(n^2) standard action. No
%   FFT/butterfly or fast-special-function algorithm is claimed.
%   [104] https://doi.org/10.1090/S0025-5718-69-99647-1
d=diagonal(:);e=offdiagonal(:);
n=numel(d);
assert(numel(e)==max(0,n-1),'Jacobi offdiagonal has incorrect length.');
if any(~isfinite(d))||any(~isfinite(e))||any(e<=0)
    error('um_gauss_basis:Jacobi','Finite diagonal and strictly positive offdiagonals required.');
end
if n==1
    Q=1;
    if nargout>1
        info=struct('nodes',d,'normalized_weights',1, ...
            'log_normalized_weights',0,'orthogonality_residual',0, ...
            'diagonalization_residual',0, ...
            'orientation','Sturm endpoint sign/positive leading coefficient');
    end
    return
end
J=diag(d)+diag(e,1)+diag(e,-1);
[V,D]=eig(J);
[nodes,perm]=sort(real(diag(D)),'ascend');
V=real(V(:,perm));
% The last component orientation sign is (-1)^(n-j) times
% the first-component orientation. Prefer whichever endpoint is
% farther from underflow/cancellation.
left=V(1,:);
right=V(end,:).*((-1).^(n-(1:n)));
anchor=left;
useLast=abs(right)>abs(left);
anchor(useLast)=right(useLast);
orientation=sign(anchor);orientation(orientation==0)=1;
V=bsxfun(@times,V,orientation);
Q=V.';
if nargout>1
    logw=um_gauss_log_weights(nodes,d,e);
    normw=exp(logw); % true Christoffel weights, possibly underflow to zero
    info=struct('nodes',nodes(:),'normalized_weights',normw(:), ...
        'log_normalized_weights',logw, ...
        'orthogonality_residual',norm(Q'*Q-eye(n),'fro'), ...
        'diagonalization_residual', ...
        norm(Q*J*Q'-diag(nodes),'fro'), ...
        'orientation','Sturm endpoint sign/positive leading coefficient');
end
end
