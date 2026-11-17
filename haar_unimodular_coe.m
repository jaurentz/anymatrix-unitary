function [U,properties] = haar_unimodular_coe(n,seed)
%HAAR_UNIMODULAR_COE Unimodular circular orthogonal (symmetric unitary).
%   Q is Haar on SU(n) and U=Q*Q.' lies in SU(n)/SO(n).
%   The transpose is non-conjugating. Do NOT replace it with Q'.
%   This is a complex symmetric unitary, not a real orthogonal matrix.
%   Ref: S. Nishigaki, PTEP 2024 (2024), 021A02.
%   https://doi.org/10.1093/ptep/ptae018
properties={'square','unitary','normal','complex','symmetric','scalable','random'};
if nargin==0,U=[];return;end
if nargin<2,seed=[];end
Q=haar_special_unitary(n,seed);
U=Q*Q.';
end
