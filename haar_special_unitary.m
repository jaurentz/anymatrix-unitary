function [U,properties] = haar_special_unitary(n,seed,output)
%HAAR_SPECIAL_UNITARY Haar-distributed SU(n) matrix (determinant one).
%   U=HAAR_SPECIAL_UNITARY(N,SEED) reproduces a sample with RNG restored.
%   U=HAAR_SPECIAL_UNITARY(N,[], 'spectral') returns the exact ENSEMBLE
%   mean eigenphase density (nonuniform for finite N>=2), not a sample.
%   A Haar U(n) draw Q is mapped to diag(conj(det(Q)),1,...,1)*Q.
%   Ref: S. Nishigaki, PTEP 2024, 021A02.
%   https://doi.org/10.1093/ptep/ptae018
%   M. Fasi and L. Robol, LAA 620 (2021), 297-321.
%   https://doi.org/10.1016/j.laa.2021.02.031
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,seed=[];end
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('haar_special_unitary',n);
    return
end
if nargin>=3
    error('haar_special_unitary:Mode','Only terminal spectral mode is supported.');
end
Q=haar_unitary(n,seed);
phase=conj(det(Q));
U=Q;
U(1,:)=phase*U(1,:);
end
