function [U,properties] = haar_unimodular_cse(n,seed)
%HAAR_UNIMODULAR_CSE Unimodular circular symplectic ensemble, order 2*n.
%   With J=[0 I;-I 0], Q Haar in SU(2*n), and the symplectic dual
%   Q^D=J*Q.'*J', form U=Q*Q^D. Then U is unitary and self-dual,
%   with determinant one and doubly degenerate eigenvalues.
%   N is the number of distinct eigenphases; matrix order is 2*N.
%   This is NOT a Haar draw from the compact symplectic group Sp(n).
%   Ref: S. Nishigaki, PTEP 2024 (2024), 021A02.
%   https://doi.org/10.1093/ptep/ptae018
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,seed=[];end
Q=haar_special_unitary(2*n,seed);
I=speye(n);
J=[sparse(n,n),I;-I,sparse(n,n)];
U=Q*J*Q.'*J';
end
