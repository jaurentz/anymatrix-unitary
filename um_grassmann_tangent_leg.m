function [T,A,Ralign,S,Delta,Q0,Q1]=um_grassmann_tangent_leg(Y,Z)
%UM_GRASSMANN_TANGENT_LEG Exact intrinsic Grassmann Levi-Civita link.
%   Y,Z are n-by-r orthonormal frames with invertible overlap.
%   Ralign=polar(Z'*Y), aligned target Zend=Z*Ralign, and
%   C=Y'*Zend=positive square root of (Z'*Y)'*(Z'*Y).
%   Principal log Delta=(Zend-Y*C)*V*diag(theta/sin(theta))*V'
%   where C=V*diag(cos(theta))*V'. Omega=Delta*Y'-Y*Delta'.
%   A=expm(Omega), A*Y=Zend, A transports every horizontal tangent
%   along the minimizing Grassmann geodesic.
%   The intrinsic tangent fiber T is
%      vec(Bnew)=kron(conj(Ralign),Q1'*A*Q0)*vec(B),
%   Q0=null(Y'), Q1=null(Z'), B of size (n-r)-by-r.
n=size(Y,1);r=size(Y,2);
[U,~,V]=svd(Z'*Y);
Ralign=U*V';
Zend=Z*Ralign;
C=Y'*Zend;
C=(C+C')/2;
[V,d]=eig(C);
lambda=real(diag(d));
if any(lambda<=1e-10)||any(lambda>1+1e-8)
    error('um_grassmann_tangent_leg:CutLocus', ...
        'Unresolved Grassmann principal angles.');
end
theta=acos(min(1,max(-1,lambda)));
scaling=ones(size(theta));
mask=theta>1e-7;
scaling(mask)=theta(mask)./sin(theta(mask));
scaling(~mask)=1+theta(~mask).^2/6;
Delta=(Zend-Y*C)*V*diag(scaling)*V';
Omega=Delta*Y'-Y*Delta';
A=expm(Omega);
Q0=null(Y');
Q1=null(Z');
S=Q1'*A*Q0;
T=kron(conj(Ralign),S);
dim=r*(n-r);
if norm(A*Y-Zend,'fro')>2e-9*n ...
    || norm(T'*T-eye(dim),'fro')>5e-9*dim
    error('um_grassmann_tangent_leg:Numerics', ...
      'Failed minimizing geodesic or tangent transport unitary test.');
end
end
