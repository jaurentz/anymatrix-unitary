function [U,properties] = jaynes_cummings(n,g,delta,t,omega_c)
%JAYNES_CUMMINGS Exact finite Jaynes-Cummings block propagator.
%   U=JAYNES_CUMMINGS(N,G,DELTA,T,OMEGA_C) acts in the tensor basis
%   |g,0>,|e,0>,|g,1>,|e,1>,...,|g,N>,|e,N>.
%   N is the MAXIMUM photon number, so order(U)=2*(N+1).
%   DELTA=omega_atom-OMEGA_C and hbar=1. Parameters are real.
%
%   Hamiltonian H=OMEGA_C*a'*a + omega_atom/2*sigma_z
%                 +G*(a*sigma_+ +a'*sigma_-).
%   Here H is COMPRESSED to photon numbers 0,...,N before exponentiation;
%   the finite propagator is exactly unitary in exact arithmetic.
%
%   Invariant internal excitation blocks are |e,k>,|g,k+1>,
%   k=0,...,N-1. In each block, the energy center is
%   omega_c*(k+1/2), with eigenenergies center +/- sqrt((delta/2)^2
%   + G^2*(k+1)). The two edge energies are -omega_atom/2 for |g,0>
%   and omega_c*N+omega_atom/2 for |e,N>.
%
%   Constructs the sparse two-state rotations analytically, including
%   zero-coupling/zero-detuning cases; no large dense expm required.
%   Defaults: G=.8, DELTA=.35, T=1, OMEGA_C=1.
%
%   Reference: Jaynes and Cummings, Proc. IEEE 51 (1963), 89-109.
properties={'square','unitary','normal','complex','sparse','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','nonnegative','finite'});
if nargin<2 || isempty(g),g=.8;end
if nargin<3 || isempty(delta),delta=.35;end
if nargin<4 || isempty(t),t=1;end
if nargin<5 || isempty(omega_c),omega_c=1;end
validateattributes(g,{'numeric'},{'scalar','real','finite'});
validateattributes(delta,{'numeric'},{'scalar','real','finite'});
validateattributes(t,{'numeric'},{'scalar','real','finite'});
validateattributes(omega_c,{'numeric'},{'scalar','real','finite'});
d=2*(n+1);
U=spalloc(d,d,4*n+2);
omega_atom=omega_c+delta;
U(1,1)=exp(1i*t*omega_atom/2);
U(d,d)=exp(-1i*t*(omega_c*n+omega_atom/2));
for k=0:n-1
    i=2*k+2;
    j=i+1;
    coupling=g*sqrt(k+1);
    split=hypot(delta/2,coupling);
    if split==0
        sinc_factor=t;
    else
        sinc_factor=sin(t*split)/split;
    end
    center=omega_c*(k+.5);
    phase=exp(-1i*t*center);
    c=cos(t*split);
    U(i,i)=phase*(c-1i*delta*sinc_factor/2);
    U(j,j)=phase*(c+1i*delta*sinc_factor/2);
    U(i,j)=phase*(-1i*coupling*sinc_factor);
    U(j,i)=U(i,j);
end
end
