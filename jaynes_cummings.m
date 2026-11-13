function [U,properties] = jaynes_cummings(n,g,delta,t,omega_c,output)
%JAYNES_CUMMINGS Exact finite Jaynes-Cummings block propagator.
%   U=JAYNES_CUMMINGS(N,G,DELTA,T,OMEGA_C) acts in the tensor basis
%   |g,0>,|e,0>,|g,1>,|e,1>,...,|g,N>,|e,N>.
%   N is the MAXIMUM photon number, order(U)=2*(N+1).
%   DELTA=omega_atom-OMEGA_C and hbar=1. Real parameters.
%
%   H=omega_c*a'*a + omega_atom/2*sigma_z
%     +G*(a*sigma_+ +a'*sigma_-).
%   The Hermitian Hamiltonian is truncated first, before exponentiation.
%   Blocks span |e,k>,|g,k+1>, k=0,...,N-1. Their energy center is
%   omega_c*(k+1/2), with splitting
%   sqrt((delta/2)^2 + G^2*(k+1)).
%   The two boundary eigenenergies are -omega_atom/2 for |g,0>
%   and omega_c*N+omega_atom/2 for the artificial cutoff |e,N>.
%
%   By default U is a sparse matrix with O(N) nonzero entries.
%   Final OUTPUT='operator' returns an O(N)-storage function handle
%   which applies independent 2x2 excitation blocks directly.
%   Both modes share exactly the same analytic block coefficients.
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
if nargin<6,output='matrix';end
implicit=um_is_operator(output);
N=2*(n+1);
omega_atom=omega_c+delta;
edge1=exp(1i*t*omega_atom/2);
edge2=exp(-1i*t*(omega_c*n+omega_atom/2));
a=complex(zeros(n,1));
b=complex(zeros(n,1));
d=complex(zeros(n,1));
for k=0:n-1
    coupling=g*sqrt(k+1);
    split=hypot(delta/2,coupling);
    if split==0
        sinc_factor=t;
    else
        sinc_factor=sin(t*split)/split;
    end
    phase=exp(-1i*t*omega_c*(k+.5));
    c=cos(t*split);
    a(k+1)=phase*(c-1i*delta*sinc_factor/2);
    b(k+1)=phase*(-1i*coupling*sinc_factor);
    d(k+1)=phase*(c+1i*delta*sinc_factor/2);
end
op=@(X,flag)um_jaynes_apply(X,flag,a,b,d,edge1,edge2);
if implicit
    U=op;
    return
end
U=spalloc(N,N,4*n+2);
U(1,1)=edge1;
U(N,N)=edge2;
for k=1:n
    i=2*k;
    U(i,i)=a(k);
    U(i,i+1)=b(k);
    U(i+1,i)=b(k);
    U(i+1,i+1)=d(k);
end
end
