function [U,properties] = vibron_2d_esqpt(N,xi,time,ell,mode)
%VIBRON_2D_ESQPT Unitary vibron model with angular-momentum sectors.
%   U=VIBRON_2D_ESQPT(N,XI,TIME) returns the sparse block-diagonal
%   unitary exp(-1i*TIME*H) for N scalar/sigma and +/- bosons (U(3)
%   two-dimensional vibron model), with all angular momenta ell=-N:N.
%   Matrix order is (N+1)*(N+2)/2. Each fixed-ell block has dimension
%   floor((N-abs(ell))/2)+1. Blocks are ordered ell=-N,...,N and
%   within each block by n_tau=abs(ell),abs(ell)+2,...,N.
%   U=VIBRON_2D_ESQPT(N,XI,TIME,ELL) selects one integer ell sector.
%
%   Hamiltonian: H=(1-XI)*n_tau+(XI/N)*(N*(N+1)*I-W^2),
%   W^2=(Dplus*Dminus+Dminus*Dplus)/2+ell^2,
%   Dplus=sqrt(2)*(tau_plus^dagger*sigma-sigma^dagger*tau_minus)
%   and Dminus=Dplus^dagger. For n_tau=t, ell fixed, n_sigma=N-t,
%   n_plus=(t+ell)/2, n_minus=(t-ell)/2:
%     (W^2)_(t,t) = 2*(N-t)*t+2*(N-t)+t+ell^2,
%     (W^2)_(t,t+2) = -2*sqrt((N-t)*(N-t-1)*
%                                     (n_plus+1)*(n_minus+1)).
%   Hence H_ell is a REAL SYMMETRIC TRIDIAGONAL matrix in this basis.
%   For full N the propagator is sparse block diagonal; within each
%   symmetry sector the propagator is generally dense.
%
%   XI in [0,1] controls the u(2)-to-so(3) dynamical symmetry change.
%   The mean-field ground-state critical value is XI=1/5. For XI>1/5
%   an excited-state quantum phase transition separates regimes of
%   near-degeneracy. Gaps between ell=0 and ell=1 levels below the ESQPT
%   typically decay algebraically in N; the corresponding 1D Lipkin
%   parity-sector gaps can decay exponentially. Neither statement
%   guarantees any degeneracy in an arbitrary finite matrix.
%
%   MODE='hamiltonian' returns the sparse real symmetric H or H_ell.
%   MODE='spectral' provides exact finite eigenphases for XI=0, XI=1,
%   and TIME=0; generic finite spectral data are marked unknown.
%   Modes 'matrix' and 'hamiltonian' preserve sector block ordering.
%   At XI=1, the so(3) Casimir labels omega=N,N-2,... >=abs(ell)
%   give E=(N*(N+1)-omega*(omega+1))/N. XI=0 gives E=n_tau.
%
%   J. Khalouf-Rivera, Q. Wang, L. F. Santos, J.-E. Garcia-Ramos,
%   M. Carvajal and F. Perez-Bernal, Degeneracy in excited-state
%   quantum phase transitions of two-level bosonic models and its
%   influence on system dynamics, Phys. Rev. A 109 (2024), 062219.
%   https://doi.org/10.1103/PhysRevA.109.062219, Eqs. (5)-(7).
properties={'square','unitary','normal','complex','scalable', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(N,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(xi),xi=.5;end
if nargin<3 || isempty(time),time=.1;end
if nargin<4,ell=[];end
if nargin<5 || isempty(mode),mode='matrix';end
validateattributes(xi,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(time,{'numeric'},{'scalar','real','finite'});
if ~isempty(ell)
    validateattributes(ell,{'numeric'}, ...
        {'scalar','integer','real','finite','>=',-N,'<=',N});
end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('vibron_2d_esqpt:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','hamiltonian','spectral'}))
    error('vibron_2d_esqpt:Mode', ...
          'Supported modes: matrix, hamiltonian, spectral.');
end
if isempty(ell)
    labels=-N:N;
    dim=(N+1)*(N+2)/2;
else
    labels=ell;
    dim=floor((N-abs(ell))/2)+1;
end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',dim);
    U.family='vibron_2d_esqpt';
    U.reference='Khalouf-Rivera et al., Phys. Rev. A 109 (2024), 062219';
    U.notes=['U(3) two-dimensional vibron model with conserved angular ' ...
        'momentum ell. A fixed sector is tridiagonal Hermitian. ' ...
        'The ESQPT near-degeneracies are finite-size spectral features, ' ...
        'not an exact spectral density formula.'];
    if time==0 || xi==0 || xi==1
        zfun=@()vibron_2d_exact_phases(N,xi,time,labels);
        U.finite.kind='exact';
        U.finite.eigenvalues=zfun;
        U.finite.moments=@(k)vibron_2d_exact_moments(zfun,k);
    end
    return
end
if strcmp(mode,'matrix')
    if isempty(ell) && N>128
        error('vibron_2d_esqpt:ExplicitOrder', ...
             'Use a selected angular-momentum sector or hamiltonian mode.');
    elseif dim>512 && ~isempty(ell)
        error('vibron_2d_esqpt:ExplicitOrder', ...
             'Sector exponential exceeds limit; use hamiltonian mode.');
    end
end
U=sparse(0,0);
for l=labels
    t=(abs(l):2:N)';
    np=(t+l)/2;
    nm=(t-l)/2;
    ns=N-t;
    diagW=2*ns.*t+2*ns+t+l^2;
    diagH=(1-xi)*t+(xi/N)*(N*(N+1)-diagW);
    offH=(2*xi/N)*sqrt(ns(1:end-1).*(ns(1:end-1)-1).* ...
                            (np(1:end-1)+1).*(nm(1:end-1)+1));
    d=numel(t);
    H=spdiags(diagH,0,d,d);
    if d>1
        H=H+sparse(1:d-1,2:d,offH,d,d);
        H=H+sparse(2:d,1:d-1,offH,d,d);
    end
    if strcmp(mode,'matrix')
        block=expm(-1i*time*full(H));
        if isempty(ell) && d<dim,block=sparse(block);end
    else
        block=H;
    end
    U=blkdiag(U,block); %#ok<AGROW>
end
end

function z=vibron_2d_exact_phases(N,xi,time,labels)
z=complex(zeros(0,1));
for l=labels
    if time==0
        energies=zeros(floor((N-abs(l))/2)+1,1);
    elseif xi==0
        energies=(abs(l):2:N)';
    else
        omega=(N:-2:abs(l))';
        energies=(N*(N+1)-omega.*(omega+1))/N;
    end
    z=[z;exp(-1i*time*energies)]; %#ok<AGROW>
end
end

function m=vibron_2d_exact_moments(zfun,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
z=zfun();
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),size(k));
end
