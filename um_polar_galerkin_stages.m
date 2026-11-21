function [g1,Ui,info]=um_polar_galerkin_stages(g0,Omega,h,stages,tol,maxit)
%UM_POLAR_GALERKIN_STAGES Source-matched constrained polar stage equations.
%   [G1,UI,INFO]=UM_POLAR_GALERKIN_STAGES(G0,OMEGA,H,STAGES,TOL,MAXIT)
%   solves only the *kinematic equations* (24b,c) of Shen, Tran,
%   Leok, J. Comput. Dyn. 9 (2022), 529--551:
%
%      Ui = P(g0 + h sum_j a_ij Uj Omega_j),
%      g1 = P(g0 + h sum_i b_i Ui Omega_i),
%
%   where P denotes the ORTHOGONAL polar factor of a matrix in
%   GL_+(n), and a_ij,b_i are s=1,2,3 Gauss--Legendre coefficients
%   as in Table 1. OMEGA_j are supplied real SKEW matrices. They
%   may be a single n-by-n skew matrix (same Ω for all stages) or
%   n-by-n-by-s. STAGES=1,2,3 (defaults 2).
%
%   IMPORTANT: THIS IS **NOT** THE FULL VPD INTEGRATOR.
%   The full symplectic variational polar algorithm ALSO solves
%   cotangent/multiplier equations (24a,d,e,f) in tandem with
%   the stages. In particular, supplied Ω_j here are not
%   guaranteed to be the Hamiltonian derivatives (24f).
%   Thus neither symplecticity nor reduced Lie--Poisson structure
%   nor any sixth-order *full dynamical method* is claimed.
%
%   On each simultaneous fixed-point iteration the source stage
%   polars are recomputed from all previous-stage values. We use
%   an SVD polar with positive orientation; if the argument leaves
%   GL_+(n) or is ill-conditioned, FAIL rather than silently switch
%   to a determinant-corrected nearest orthogonal matrix.
%   The max stage fixed-point residual is reevaluated after
%   convergence. Max iterations default 200, tolerance 1e-12.
%
%   The returned G1 and every UI(:,:,j) lie in SO(n) up to
%   floating-point arithmetic. They are rotations, not a
%   symplectic phase-space Jacobian. A left change of spatial
%   frame S in SO(n) gives G1new=S*G1 and Uinew=S*Ui,
%   since P(S*A)=S*P(A). Each polar uses dense SVD O(n^3),
%   each full iteration O(stages^2*n^3) for generic matrices.
%
%   Source: Shen--Tran--Leok (2022), Eqs. (24b,c), Table 1.
%   https://doi.org/10.3934/jcd.2022003
if ~(isnumeric(g0)&&isreal(g0)&&ismatrix(g0) ...
        &&size(g0,1)==size(g0,2)&&size(g0,1)>=2 ...
        &&size(g0,1)<=32&&all(isfinite(g0(:))))
    error('um_polar_galerkin_stages:Configuration', ...
        'G0 must be a finite real n-by-n matrix, 2<=n<=32.');
end
n=size(g0,1);
if norm(g0'*g0-eye(n),'fro')>2e-10*n || det(g0)<0
    error('um_polar_galerkin_stages:Configuration', ...
        'G0 must already lie in SO(n); no silent projection.');
end
if nargin<4||isempty(stages),stages=2;end
[A,b,c,order]=um_polar_galerkin_coefficients(stages);
if nargin<5||isempty(tol),tol=1e-12;end
if nargin<6||isempty(maxit),maxit=200;end
if ~(isnumeric(tol)&&isscalar(tol)&&isreal(tol) ...
        &&isfinite(tol)&&tol>=1e-15&&tol<=1e-4)
    error('um_polar_galerkin_stages:Tolerance', ...
        'TOL must lie in [1e-15,1e-4].');
end
if ~(isnumeric(maxit)&&isscalar(maxit)&&isreal(maxit) ...
        &&isfinite(maxit)&&maxit==fix(maxit) ...
        &&maxit>=1&&maxit<=10000)
    error('um_polar_galerkin_stages:Iterations', ...
        'MAXIT must be integer 1..10000.');
end
if ~(isnumeric(h)&&isscalar(h)&&isreal(h)&&isfinite(h))
    error('um_polar_galerkin_stages:Step','h must be real finite.');
end
if ~(isnumeric(Omega)&&isreal(Omega)&& ...
        size(Omega,1)==n&&size(Omega,2)==n ...
        &&ndims(Omega)<=3&&all(isfinite(Omega(:))))
    error('um_polar_galerkin_stages:Omega', ...
        'OMEGA must be real n-by-n or n-by-n-by-stages.');
end
if size(Omega,3)==1
    Omega=repmat(Omega,1,1,stages);
elseif size(Omega,3)~=stages
    error('um_polar_galerkin_stages:Omega', ...
        'Expected one or STAGES skew matrices.');
end
for j=1:stages
    W=Omega(:,:,j);
    if norm(W+W','fro')>2e-12*max(1,norm(W,'fro'))
        error('um_polar_galerkin_stages:Omega', ...
            'Every Ω_j must be real skew-symmetric.');
    end
end
Ui=repmat(g0,1,1,stages);
residual=Inf;
iterations=0;
for k=1:maxit
    next=zeros(n,n,stages);
    for i=1:stages
        Y=g0;
        for j=1:stages
            Y=Y+h*A(i,j)*Ui(:,:,j)*Omega(:,:,j);
        end
        next(:,:,i)=source_polar(Y);
    end
    iterateDifference=0;
    for j=1:stages
        iterateDifference=max(iterateDifference, ...
            norm(next(:,:,j)-Ui(:,:,j),'fro'));
    end
    Ui=next;
    iterations=k;
    if iterateDifference<tol
        break
    end
end
% Recompute the actual source fixed-point residual from final stages.
residual=0;
for i=1:stages
    Y=g0;
    for j=1:stages
        Y=Y+h*A(i,j)*Ui(:,:,j)*Omega(:,:,j);
    end
    residual=max(residual,norm(Ui(:,:,i)-source_polar(Y),'fro'));
end
if residual>=tol
    error('um_polar_galerkin_stages:NoConvergence', ...
        'Polar stage fixed-point residual %.3e exceeds TOL %.3e.', ...
        residual,tol);
end
Y=g0;
for j=1:stages
    Y=Y+h*b(j)*Ui(:,:,j)*Omega(:,:,j);
end
g1=source_polar(Y);
maxOrth=max(norm(g1'*g1-eye(n),'fro'),0);
for j=1:stages
    maxOrth=max(maxOrth,norm(Ui(:,:,j)'*Ui(:,:,j)-eye(n),'fro'));
end
info=struct('A',A,'b',b,'c',c,'gauss_order',order, ...
    'iterations',iterations,'stage_residual',residual, ...
    'iteration_difference',iterateDifference, ...
    'maximum_orthogonality_defect',maxOrth, ...
    'endpoint_orthogonality_defect',norm(g1'*g1-eye(n),'fro'), ...
    'is_full_variational_integrator',false, ...
    'has_momentum_or_symplectic_update',false, ...
    'source_equations','Shen--Tran--Leok (2022), Eqs. (24b,c)', ...
    'notes',['Source-exact polar-projected stage and endpoint ', ...
        'kinematics for prescribed skew Ω. Full symplectic VPD ', ...
        'requires the coupled (24a,d,e,f) equations, not solved here.']);
end

function Q=source_polar(A)
n=size(A,1);
if any(~isfinite(A(:)))||det(A)<=0
    error('um_polar_galerkin_stages:PolarDomain', ...
        'Polar argument must belong to GL_+(n).');
end
[U,S,V]=svd(A);
sigma=diag(S);
if min(sigma)<=1e-12*max(sigma)
    error('um_polar_galerkin_stages:PolarDomain', ...
        'Polar argument is singular or ill-conditioned.');
end
Q=U*V';
if det(Q)<=0
    error('um_polar_galerkin_stages:PolarDomain', ...
        'Polar factor left SO(n) despite positive determinant input.');
end
end
