function [U,properties] = kicked_top(n,kappa,p)
%KICKED_TOP Quantum kicked top on the spin-j Hilbert space, n=2*j+1.
%   U = exp(-1i*kappa*Jz^2/(2*j))*exp(-1i*p*Jy), j=(n-1)/2.
%   Basis: magnetic quantum numbers m=j,j-1,...,-j.
%   n>=2 permits both integer and half-integer spin.
%   At p=0 the eigenvalues are exact diagonal twist phases.
%   At kappa=0 the eigenvalues are exp(-1i*p*m), m=-j,...,j.
%   References: F. Haake, Quantum Signatures of Chaos (2010);
%   V. M. Bastidas et al., Phys. Rev. Lett. 112 (2014), 140408.
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(kappa),kappa=3;end
if nargin<3 || isempty(p),p=pi/2;end
validateattributes(kappa,{'numeric'},{'scalar','real','finite'});
validateattributes(p,{'numeric'},{'scalar','real','finite'});
j=(n-1)/2;
m=(j:-1:-j)';
raising=diag(sqrt((1:n-1).*(n-(1:n-1))),1);
Jy=(raising-raising')/(2i);
twist=diag(exp(-1i*kappa*(m.^2)/(2*j)));
U=twist*expm(-1i*p*Jy);
end
