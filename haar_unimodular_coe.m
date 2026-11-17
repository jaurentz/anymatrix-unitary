function [U,properties] = haar_unimodular_coe(n,seed,output)
%HAAR_UNIMODULAR_COE Unimodular COE: symmetric unitary SU(n)/SO(n).
%   Haar Q in SU(n), U=Q*Q.'; transpose here is NON-conjugating.
%   U=HAAR_UNIMODULAR_COE(N,[], 'spectral') returns exact finite-N
%   ensemble density, including the exceptional N=2 |sin(theta)|/4.
%   Not Haar in the real orthogonal group.
%   Ref: S. Nishigaki, PTEP 2024, 021A02.
%   https://doi.org/10.1093/ptep/ptae018
properties={'square','unitary','normal','complex','symmetric','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,seed=[];end
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('haar_unimodular_coe',n);
    return
end
if nargin>=3
    error('haar_unimodular_coe:Mode','Only terminal spectral mode is supported.');
end
Q=haar_special_unitary(n,seed);
U=Q*Q.';
end
