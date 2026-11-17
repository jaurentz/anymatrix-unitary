function [U,properties] = haar_special_unitary(n,seed)
%HAAR_SPECIAL_UNITARY Haar-distributed SU(n) matrix (determinant one).
%   U=HAAR_SPECIAL_UNITARY(N,SEED) reproduces a draw without advancing the
%   caller's random stream. N=1 returns 1. A Haar U(n) draw Q is mapped
%   to diag(conj(det(Q)),1,...,1)*Q; the pushforward is Haar on SU(n).
%   The normalized expected eigenphase density for n>=2 is
%   [1-(-1)^n*(2/n)*cos(n*theta)]/(2*pi).
%   Ref: S. Nishigaki, PTEP 2024 (2024), 021A02.
%   https://doi.org/10.1093/ptep/ptae018
%   M. Fasi and L. Robol, sampling Haar eigenvalues, LAA (2021).
%   https://doi.org/10.1016/j.laa.2021.02.031
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
if nargin<2,seed=[];end
if isempty(seed)
    Q=haar_unitary(n);
else
    Q=haar_unitary(n,seed);
end
phase=conj(det(Q));
U=Q;
U(1,:)=phase*U(1,:);
end
