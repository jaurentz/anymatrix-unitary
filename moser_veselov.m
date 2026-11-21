function [F,properties]=moser_veselov(n,h,J,M,steps,mode)
%MOSER_VESELOV Implicit discrete generalized rigid-body SO(n) increments.
%   F=MOSER_VESELOV(N,H,J,M) returns the near-identity SO(N) body
%   increment obtained from the Moser--Veselov relation
%       H*M = F'*J - J*F,
%   with J=J' positive definite and M=-M' real body angular momentum.
%   The discrete Euler--Arnold update in this SAME convention is
%       M_next = F*M*F',   Q_next = Q*F'.
%   In particular, Q*M*Q' is exactly conserved up to rounding.
%   Note carefully the transpose of F in the orientation update;
%   substituting Q_next=Q*F while retaining the above momentum
%   equation is an inconsistent choice of frame convention.
%
%   Unlike a generic polar orthogonal factor, F is determined by
%   solving the nonlinear mechanical implicit matrix equation.
%   A damped analytic Newton method in skew Cayley coordinates selects
%   the SO(N) branch CONTINUOUS FROM IDENTITY at small H*M. Other
%   admissible distant branches may exist; arbitrary inputs need not
%   admit a solution in this local chart. Failure throws an error, not
%   a silently approximate orthogonal factor. N=3,...,16.
%
%   Defaults H=.15, deterministic positive tridiagonal J and skew
%   banded M; STEPS=1, MODE='matrix'. The physical step H is absorbed
%   into the discrete momentum in the relation above.
%
%   MODE='matrix'     First n-by-n SO(n) increment F.
%   MODE='model'      First F, J,M, next momentum and orientation,
%                     residual and explicitly defined Casimirs.
%   MODE='trajectory' STEPS consecutive implicit increments with
%                     Q(:,:,k), M(:,:,k), F(:,:,k), residuals.
%                     Q(:,:,1)=I, M(:,:,1)=initial M.
%   MODE='spectral'  Exact identity-phase law for M=0; otherwise
%                     unknown finite phase law with SO(n) pairing.
%   No matrix-free acceleration is claimed for a dense nonlinear
%   solve. SO(n) orthogonality is NOT symplecticity of a full
%   rigid-body phase-space Jacobian.
%
%   Moser and Veselov (1991), Commun. Math. Phys. 139, 217-243.
%   [138] https://doi.org/10.1007/BF02352494
%   Cardoso and Miraldo (2022), J. Comput. Appl. Math. 402, 113814.
%   [139] https://doi.org/10.1016/j.cam.2021.113814
%   Shen, Tran and Leok (2022), J. Comput. Dyn. 9, 529-551:
%   different high-order polar variational integrator, NOT coded here.
%   [140] https://doi.org/10.3934/jcd.2022003
properties={'square','unitary','orthogonal','normal','real', ...
            'parameter-dependent'};
if nargin==0,F=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',3,'<=',16});
if nargin<2 || isempty(h),h=.15;end
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
if nargin<3 || isempty(J)
    d=(1:n)';
    J=diag(1+.2*d+.08*sin(d));
    if n>1
        J=J+.04*diag(ones(n-1,1),1)+.04*diag(ones(n-1,1),-1);
    end
end
if nargin<4 || isempty(M)
    M=zeros(n);
    for k=1:n-1,M(k,k+1)=.17*cos(k);end
    for k=1:n-2,M(k,k+2)=.1*sin(k);end
    M=M-M';
end
if nargin<5 || isempty(steps),steps=1;end
if nargin<6 || isempty(mode),mode='matrix';end
validateattributes(steps,{'numeric'},{'scalar','integer','finite','>=',1,'<=',40});
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('moser_veselov:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','trajectory','spectral'}))
    error('moser_veselov:Mode','Unknown output mode.');
end
if ~isnumeric(J) || ~ismatrix(J) || ~isequal(size(J),[n n]) ...
        || ~isreal(J) || ~all(isfinite(J(:)))
    error('moser_veselov:Inertia','J must be a finite real n-by-n inertia matrix.');
end
if norm(J-J','fro')>1e-12*max(1,norm(J,'fro'))
    error('moser_veselov:Inertia','J must be symmetric.');
end
J=(J+J')/2;
[~,p]=chol(J);
if p~=0 || rcond(J)<1e-12
    error('moser_veselov:Inertia', ...
        'J must be well-conditioned, symmetric and positive definite.');
end
if ~isnumeric(M) || ~ismatrix(M) || ~isequal(size(M),[n n]) ...
        || ~isreal(M) || ~all(isfinite(M(:)))
    error('moser_veselov:Momentum','M must be a finite real n-by-n skew matrix.');
end
if norm(M+M','fro')>1e-12*max(1,norm(M,'fro'))
    error('moser_veselov:Momentum','M must be skew-symmetric.');
end
M=(M-M')/2;
if strcmp(mode,'spectral')
    empty=struct('kind','unknown','density',[],'support',[], ...
                 'atoms',[],'moments',[]);
    finite=struct('kind','unknown','eigenvalues',[], ...
                  'moments',[],'atoms',[]);
    notes=['Generic implicit SO(n) increment has no closed analytic ', ...
      'trace eigenphase law. Conjugate eigenphase pairing; det(F)=1; ', ...
      'odd n guarantees one eigenvalue +1. The unitary F spectrum ', ...
      'must not be confused with spectra of the full symplectic ', ...
      'rigid-body map.'];
    if ~any(M(:))
        finite.kind='exact';
        finite.eigenvalues=@()ones(n,1);
        finite.moments=@(k)ones(size(k));
        finite.atoms=struct('angles',0,'weights',1);
        notes=['Near-identity branch with zero body angular momentum ', ...
            'has F=I. Exact sample spectrum is an n-fold atom at 1.'];
    end
    F=struct('family','moser_veselov','dimension',n, ...
       'normalization','probability (density with respect to dtheta)', ...
       'finite',finite,'limit',empty,'ensemble',empty,'underlying',empty, ...
       'reference','Moser and Veselov (1991); Cardoso and Miraldo (2022)', ...
       'notes',notes);
    return
end
if strcmp(mode,'trajectory')
    F_steps=zeros(n,n,steps);
    Q_steps=zeros(n,n,steps+1);
    M_steps=zeros(n,n,steps+1);
    residuals=zeros(steps,1);
    Q_steps(:,:,1)=eye(n);
    M_steps(:,:,1)=M;
    for k=1:steps
        B=um_moser_veselov_solve(h,J,M_steps(:,:,k));
        F_steps(:,:,k)=B;
        M_steps(:,:,k+1)=B*M_steps(:,:,k)*B';
        Q_steps(:,:,k+1)=Q_steps(:,:,k)*B';
        residuals(k)=norm(B'*J-J*B-h*M_steps(:,:,k),'fro');
    end
    F=struct('F',F_steps,'Q',Q_steps,'M',M_steps, ...
       'J',J,'h',h,'steps',steps,'residuals',residuals, ...
       'initial_spatial_momentum',M, ...
       'notes',['Q(k+1)=Q(k)*F(k)''; M(k+1)=F(k)*M(k)*F(k)''. ', ...
                'Q(k)*M(k)*Q(k)'' is the conserved spatial momentum. ', ...
                'This does not assert that the phase-space Jacobian ', ...
                'is orthogonal.']);
    return
end
F=um_moser_veselov_solve(h,J,M);
if strcmp(mode,'model')
    Q_next=F';
    M_next=F*M*F';
    F=struct('F',F,'Q_next',Q_next,'M_next',M_next, ...
        'J',J,'M',M,'h',h, ...
        'residual',norm(F'*J-J*F-h*M,'fro'), ...
        'spatial_momentum_residual',norm(Q_next*M_next*Q_next'-M,'fro'), ...
        'determinant',det(F), ...
        'notes','SO(n) body increment; not a full orthogonal phase-space map.');
end
end
