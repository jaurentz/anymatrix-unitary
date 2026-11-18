function [U,properties,ensemble] = lohe_matrix_flow(d,nosc,kappa,t,H,U0,which,steps,mode)
%LOHE_MATRIX_FLOW Nonlinear coupled unitary matrices (Lohe/Kuramoto model).
%   U=LOHE_MATRIX_FLOW(D,NOSC,KAPPA,T,H,U0,WHICH,STEPS) returns the
%   WHICH-th D-by-D unitary matrix of an NOSC-oscillator Lohe ensemble.
%   NOSC defaults to 3, KAPPA to 1, T to 1, WHICH to 1, and STEPS to
%   max(8,ceil(32*abs(T))). H and U0 are D-by-D-by-NOSC arrays.
%   If unspecified, deterministic Hermitian H and independently seeded
%   Haar-unitary initial matrices U0 are used (caller RNG is unchanged).
%
%   The continuous-time system is
%      dU_j/dt = A_j(U) U_j,
%      A_j = -i H_j + (KAPPA/2)*(M*U_j' - U_j*M'),
%      M = (1/NOSC)*sum_l U_l.
%   H_j must be Hermitian and U_j(0) must be unitary. This is a nonlinear
%   dynamical system on U(D)^NOSC: U_j(t) is unitary for all real T.
%
%   The exponential midpoint method uses a synchronous half-step
%   predictor and a full update exp(h*A_j(Uhalf))*U_j. Each step is
%   unitary up to roundoff and second-order accurate for fixed T as
%   STEPS increases. Its trajectory is an APPROXIMATION, not an exact
%   solution of the continuous-time Lohe equations.
%
%   [U,PROPERTIES,ALL]=LOHE_MATRIX_FLOW(...) also returns the complete
%   D-by-D-by-NOSC evolved ensemble. Terminal MODE='ensemble' returns
%   that ensemble as the first output; MODE='spectral' returns an
%   unknown-law spectral descriptor (no general analytic phase density).
%   The default MODE='matrix' returns only oscillator WHICH.
%
%   Special cases: KAPPA=0 is decoupled exp(-i*T*H_j)*U_j(0);
%   NOSC=1 is exact for any number of steps; D=1 is Kuramoto (phase
%   sign is -H); identical H and synchronized U0 remain synchronized.
%   Starting in SU(D) does NOT generally guarantee det(U_j(t))=1,
%   because the Lohe coupling need not have zero trace.
%
%   M. A. Lohe, Non-Abelian Kuramoto models and synchronization,
%   J. Phys. A 42 (2009), 395101.
%   https://doi.org/10.1088/1751-8113/42/39/395101
%   H. Choi, S.-Y. Ha and H. Park, Emergent behaviors of discrete Lohe aggregation
%   flows, Discrete Contin. Dyn. Syst. B (2022).
%   https://doi.org/10.3934/dcdsb.2021308
properties={'square','unitary','normal','complex','scalable', ...
            'parameter-dependent'};
ensemble=[];
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if d>128
    error('lohe_matrix_flow:Order','Explicit order is capped at 128.');
end
if nargin<2||isempty(nosc),nosc=3;end
if nargin<3||isempty(kappa),kappa=1;end
if nargin<4||isempty(t),t=1;end
validateattributes(nosc,{'numeric'},{'scalar','integer','positive','finite','<=',64});
validateattributes(kappa,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<5,H=[];end
if nargin<6,U0=[];end
if nargin<7||isempty(which),which=1;end
if nargin<8||isempty(steps),steps=max(8,ceil(32*abs(t)));end
if nargin<9||isempty(mode),mode='matrix';end
validateattributes(which,{'numeric'},{'scalar','integer','positive','finite','<=',nosc});
validateattributes(steps,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',100000});
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('lohe_matrix_flow:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','ensemble','spectral'}))
    error('lohe_matrix_flow:Mode', ...
          'Supported modes are matrix, ensemble, spectral.');
end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('lohe_matrix_flow',d);
    U.reference='Lohe, J. Phys. A 42 (2009), 395101';
    U.notes=['Nonlinear unitary matrix flow with time-dependent discrete ' ...
        'eigenphases; no general analytic finite trace spectrum or ' ...
        'ensemble expectation is asserted. The integrator is an ' ...
        'exponential midpoint approximation that preserves unitarity.'];
    return
end
if isempty(H)
    H=zeros(d,d,nosc);
    diagvalues=linspace(-1,1,d)';
    band=diag(ones(d-1,1),1)+diag(ones(d-1,1),-1);
    D=diag(diagvalues);
    for j=1:nosc
        H(:,:,j)=(1+0.2*(j-1))*D+(0.08+0.04*j)*band ...
            +0.3*(j-(nosc+1)/2)*eye(d);
    end
else
    check_array(H,d,nosc,'H');
    for j=1:nosc
        X=H(:,:,j);
        if norm(X-X','fro')>1e-11*max(1,norm(X,'fro'))
            error('lohe_matrix_flow:Hermitian', ...
                  'Every H(:,:,j) must be Hermitian.');
        end
    end
end
if isempty(U0)
    U0=complex(zeros(d,d,nosc));
    for j=1:nosc
        U0(:,:,j)=haar_unitary(d,1900+j);
    end
else
    check_array(U0,d,nosc,'U0');
    for j=1:nosc
        X=U0(:,:,j);
        if norm(X'*X-eye(d),'fro')>1e-10*max(1,sqrt(d))
            error('lohe_matrix_flow:Unitary', ...
                  'Every U0(:,:,j) must be unitary.');
        end
    end
end
ensemble=U0;
if t~=0
    h=t/steps;
    for iteration=1:steps
        A=lohe_generators(ensemble,H,kappa,nosc);
        halfway=complex(zeros(d,d,nosc));
        for j=1:nosc
            halfway(:,:,j)=expm((h/2)*A(:,:,j))*ensemble(:,:,j);
        end
        Amid=lohe_generators(halfway,H,kappa,nosc);
        next=complex(zeros(d,d,nosc));
        for j=1:nosc
            next(:,:,j)=expm(h*Amid(:,:,j))*ensemble(:,:,j);
        end
        ensemble=next;
    end
end
if strcmp(mode,'ensemble')
    U=ensemble;
else
    U=ensemble(:,:,which);
end
end

function check_array(A,d,nosc,name)
if ~isnumeric(A)||~isfloat(A)||~isequal(size(A,1),d) ...
        ||~isequal(size(A,2),d)||~isequal(size(A,3),nosc) ...
        ||~all(isfinite(A(:)))
    error('lohe_matrix_flow:Input', ...
          '%s must be a finite D-by-D-by-NOSC numeric array.',name);
end
end

function A=lohe_generators(U,H,kappa,nosc)
[d,~,~]=size(U);
M=sum(U,3)/nosc;
A=complex(zeros(d,d,nosc));
for j=1:nosc
    X=U(:,:,j);
    Aj=-1i*H(:,:,j)+(kappa/2)*(M*X'-X*M');
    A(:,:,j)=(Aj-Aj')/2;  % removes only floating-point skew drift
end
end
