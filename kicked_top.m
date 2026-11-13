function [U,properties] = kicked_top(n,kappa,p,output)
%KICKED_TOP Quantum kicked top on the spin-j Hilbert space, n=2*j+1.
%   U = exp(-1i*kappa*Jz^2/(2*j))*exp(-1i*p*Jy), j=(n-1)/2.
%   Basis: magnetic quantum numbers m=j,j-1,...,-j.
%   n>=2 permits both integer and half-integer spin.
%   At p=0 the eigenvalues are exact diagonal twist phases.
%   At kappa=0 the eigenvalues are exp(-1i*p*m), m=-j,...,j.
%   Final OUTPUT='operator' reuses UM_KICK_TRIDIAG_APPLY for a
%   sparse-tridiagonal EXPMV spin rotation (Chebyshev fallback).
%   References: F. Haake, Quantum Signatures of Chaos (2010);
%   V. M. Bastidas et al., Phys. Rev. Lett. 112 (2014), 140408.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [28] https://doi.org/10.1007/978-3-642-05428-0
%   [29] https://doi.org/10.1103/PhysRevLett.112.140408
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(kappa),kappa=3;end
if nargin<3 || isempty(p),p=pi/2;end
validateattributes(kappa,{'numeric'},{'scalar','real','finite'});
validateattributes(p,{'numeric'},{'scalar','real','finite'});
if nargin<4,output='matrix';end
implicit=um_is_operator(output);
j=(n-1)/2;
m=(j:-1:-j)';
upper=sqrt((1:n-1).*(n-(1:n-1)))/(2i);
phase=exp(-1i*kappa*(m.^2)/(2*j));
if implicit
    U=@(X,flag)um_kick_tridiag_apply(X,flag,phase,zeros(n,1),upper(:),p);
    return
end
% The Hermitian spin generator is sparse; EXPM produces a dense rotation.
Jy=um_sparse_hermitian_tridiag(zeros(n,1),upper(:));
U=bsxfun(@times,phase,full(expm(-1i*p*Jy)));
end
