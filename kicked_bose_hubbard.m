function [U,properties] = kicked_bose_hubbard(n,J,g,tau_hop,tau_kick,bias,output)
%KICKED_BOSE_HUBBARD Floquet operator of a two-mode Bose-Hubbard dimer.
%   U=KICKED_BOSE_HUBBARD(N,J,G,TAU_HOP,TAU_KICK,BIAS)
%   uses N identical bosons in two sites, in the fixed-number Fock basis
%   |k,N-k>, k=0,...,N. The resulting unitary has order N+1, NOT 2^N.
%
%   U=exp(-i*TAU_KICK*Hdiag)*exp(-i*TAU_HOP*Hhop), where
%       Hhop=-J*(a1'*a2+a2'*a1),
%       Hdiag=G/2*(n1*(n1-1)+n2*(n2-1))+BIAS/2*(n1-n2).
%   The first factor is diagonal and the second is a unitary exponential
%   of a real symmetric tridiagonal matrix. These Floquet factors are
%   evaluated explicitly, not by truncating a nonunitary propagator.
%   With final OUTPUT='operator', return a matrix-free function handle
%   for both forward and conjugate-adjoint actions. The hopping action
%   uses the shared tridiagonal Chebyshev-Bessel exponential kernel.
%
%   Defaults: J=1, G=.4, TAU_HOP=1, TAU_KICK=1, BIAS=.17.
%   With G=BIAS=0 the exact eigenvalues are
%       exp(-i*TAU_HOP*J*(2*k-N)), k=0,...,N.
%   With J=0 the exact eigenvalues are the diagonal kick phases.
%   N=0 is the one-dimensional vacuum sector.
%
%   References: Strzys, Graefe, Korsch, New J. Phys. 10 (2008), 013024;
%   Khripkov, Cohen, Vardi, Phys. Rev. E 87 (2013), 012910.
properties={'square','unitary','normal','complex','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','nonnegative','finite'});
if nargin<2 || isempty(J),J=1;end
if nargin<3 || isempty(g),g=.4;end
if nargin<4 || isempty(tau_hop),tau_hop=1;end
if nargin<5 || isempty(tau_kick),tau_kick=1;end
if nargin<6 || isempty(bias),bias=.17;end
validateattributes(J,{'numeric'},{'scalar','real','finite'});
validateattributes(g,{'numeric'},{'scalar','real','finite'});
validateattributes(tau_hop,{'numeric'},{'scalar','real','finite'});
validateattributes(tau_kick,{'numeric'},{'scalar','real','finite'});
validateattributes(bias,{'numeric'},{'scalar','real','finite'});
if nargin<7,output='matrix';end
implicit=um_is_operator(output);
k=(0:n)';
E=.5*g*(k.*(k-1)+(n-k).*(n-k-1))+.5*bias*(2*k-n);
if n==0
    if implicit
        U=@(X,flag)um_kick_tridiag_apply(X,flag,1,0,[],tau_hop);
    else
        U=1;
    end
    return
end
hop=sqrt(((1:n).*(n-(1:n)+1))).';
if implicit
    phase=exp(-1i*tau_kick*E);
    upper=-J*hop;
    U=@(X,flag)um_kick_tridiag_apply(X,flag,phase,zeros(n+1,1),upper,tau_hop);
    return
end
Hhop=-J*(diag(hop,1)+diag(hop,-1));
U=bsxfun(@times,exp(-1i*tau_kick*E),expm(-1i*tau_hop*Hhop));
end
