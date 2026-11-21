function C=um_wigner_real_basis(ell)
%UM_WIGNER_REAL_BASIS Unitary Condon-Shortley complex-to-real SH basis.
%   C maps real-harmonic expansion coefficients into complex Y_l^m
%   coefficients, using m=-ell:ell in BOTH coefficient orderings.
%   m>0 labels cosine harmonics; m<0 labels sine harmonics:
%      R_l^m  = (Y_l^(-m)+(-1)^m Y_l^m)/sqrt(2), m>0;
%      R_l^(-m)= (Y_l^(-m)-(-1)^m Y_l^m)/(i*sqrt(2)).
%   These are real functions under Y_l^(-m)=(-1)^m conj(Y_l^m).
%   The change of basis C has exactly 4*ell+1 nonzeros and is sparse.
%   The rotation in the real harmonic basis is C'*D_complex*C,
%   real orthogonal in exact arithmetic. Naively REAL(D_complex)
%   does NOT give that matrix.
validateattributes(ell,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
n=2*ell+1;
C=spalloc(n,n,4*ell+1);
C(ell+1,ell+1)=1;
for m=1:ell
    neg=ell+1-m;pos=ell+1+m;
    C(neg,pos)=1/sqrt(2);
    C(pos,pos)=(-1)^m/sqrt(2);
    C(neg,neg)=-1i/sqrt(2);
    C(pos,neg)=1i*(-1)^m/sqrt(2);
end
end
