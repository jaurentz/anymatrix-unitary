function [U,properties] = extended_lipkin_propagator(n,lambda,alpha,t,mode)
%EXTENDED_LIPKIN_PROPAGATOR Collective two-level nuclear-boson evolution.
%   U=EXTENDED_LIPKIN_PROPAGATOR(N,LAMBDA,ALPHA,T) forms the (N+1)-by-(N+1)
%   unitary U=expm(-1i*T*H) for the fixed-N two-level boson sector.
%   With n_t=t_boson^dagger*t_boson and Q(alpha)=s^dagger*t+t^dagger*s
%   +alpha*n_t, use H=(1-lambda)*n_t-(lambda/N)*Q(alpha)^2.
%   Note that 't' is physical evolution time, not the t-boson operator.
%   The basis is ordered by occupation of the t boson: k=0,...,N.
%   In this basis n_t=diag(k), Q has diagonal alpha*k and real
%   offdiagonals sqrt((k+1)*(N-k)), so H is real symmetric
%   pentadiagonal. The propagator is generally dense and unitary.
%
%   LAMBDA in [0,1], ALPHA real (default 0), T real (default 1).
%   The mean-field critical line is lambda=1/(5+alpha^2).
%   At ALPHA=0 the standard Lipkin model has a second-order transition;
%   for nonzero ALPHA the transition is generically first order.
%   These phase-transition claims concern the thermodynamic limit.
%
%   MODE='hamiltonian' returns the sparse symmetric H rather than U.
%   MODE='spectral' gives exact finite eigenphases for LAMBDA=0,
%   for LAMBDA=1 and ALPHA=0, and for T=0; other parameters are
%   represented as unknown, without performing a dense eigensolve.
%   The default matrix mode is intended for moderate N; for large N
%   use 'hamiltonian' with a structured matrix-function method.
%
%   Reference: S. Baid, A. Saiz, L. Lamata, P. Perez-Fernandez,
%   A. M. Romero, A. Rios, J. M. Arias and J. E. Garcia-Ramos,
%   Extended Lipkin model: Proposal for implementation in a quantum
%   platform and machine learning analysis of its phase diagram,
%   Phys. Rev. C 110 (2024), 044318, Eqs. (4)-(5).
%   https://doi.org/10.1103/PhysRevC.110.044318
properties={'square','unitary','normal','complex','scalable', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(lambda),lambda=.2;end
if nargin<3 || isempty(alpha),alpha=0;end
if nargin<4 || isempty(t),t=1;end
validateattributes(lambda,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<5 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('extended_lipkin_propagator:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','hamiltonian','spectral'}))
    error('extended_lipkin_propagator:Mode', ...
          'Supported modes are matrix, hamiltonian, spectral.');
end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('kicked_bose_hubbard',n);
    U.family='extended_lipkin_propagator';
    U.reference='Baid et al., Phys. Rev. C 110 (2024), 044318';
    U.notes=['Unitary evolution of a fixed-N extended Lipkin Hamiltonian. ' ...
        'Mean-field critical lambda=1/(5+alpha^2) pertains to the ' ...
        'thermodynamic limit; finite matrix has discrete spectrum.'];
    zfun=[];
    if t==0
        zfun=@()ones(n+1,1);
    elseif lambda==0
        zfun=@()exp(-1i*t*(0:n)');
    elseif lambda==1 && alpha==0
        zfun=@()exp((4i*t/n)*((0:n)'-n/2).^2);
    end
    if ~isempty(zfun)
        U.finite.kind='exact';
        U.finite.eigenvalues=zfun;
        U.finite.moments=@(k)elp_finite_moments(zfun,k);
    end
    return
end
k=(0:n)';
b=sqrt((1:n)'.*(n-(1:n)'+1));
ix=(1:n)';
Q=sparse([ix;ix+1],[ix+1;ix],[b;b],n+1,n+1);
Q=Q+spdiags(alpha*k,0,n+1,n+1);
H=(1-lambda)*spdiags(k,0,n+1,n+1)-(lambda/n)*(Q*Q);
if strcmp(mode,'hamiltonian')
    U=H;
    return
end
if n>512
    error('extended_lipkin_propagator:ExplicitOrder', ...
          'Use hamiltonian mode for large N; explicit exponential is dense.');
end
U=expm(-1i*t*full(H));
end

function m=elp_finite_moments(zfun,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
z=zfun();
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),size(k));
end
