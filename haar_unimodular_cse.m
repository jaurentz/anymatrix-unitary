function [U,properties] = haar_unimodular_cse(n,seed,output)
%HAAR_UNIMODULAR_CSE Self-dual unimodular CSE, matrix order 2*n.
%   Q is Haar in SU(2n), J=[0 I;-I 0], U=Q*J*Q.'*J'.
%   Eigenvalues have exact double degeneracy (Kramers pairing).
%   U=HAAR_UNIMODULAR_CSE(N,[], 'spectral') returns finite-N mean density.
%   For N=1 the matrix is I_2 and the density is atomic, not continuous.
%   Ref: S. Nishigaki, PTEP 2024, 021A02.
%   https://doi.org/10.1093/ptep/ptae018
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,seed=[];end
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('haar_unimodular_cse',n);
    return
end
if nargin>=3
    error('haar_unimodular_cse:Mode','Only terminal spectral mode is supported.');
end
Q=haar_special_unitary(2*n,seed);
I=speye(n);
J=[sparse(n,n),I;-I,sparse(n,n)];
U=Q*J*Q.'*J';
end
