function [U,properties] = polar_unitary(n,A,output)
%POLAR_UNITARY Unitary polar factor of a nonsingular N-by-N matrix.
%   U=POLAR_UNITARY(N,A) returns L*R' from the singular value
%   decomposition A=L*S*R'. A must be nonsingular, so the unitary
%   polar factor is uniquely determined. For real A, U is orthogonal.
%   Default A=I+0.4*upper_shift is nonsingular and nonnormal.
%   Reference: N. J. Higham, SIAM J. Sci. Stat. Comput. 7 (1986),
%   1160-1174. https://doi.org/10.1137/0907079
properties={'square','unitary','normal','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3
    if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
        U=unitary_spectral_reference('polar_unitary',n,A);
        return
    end
    error('polar_unitary:Mode','Only terminal spectral mode is supported for this generator.');
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(A)
    A=eye(n)+0.4*diag(ones(n-1,1),1);
end
validateattributes(A,{'numeric'},{'size',[n n],'finite'});
[L,S,R]=svd(full(A));
s=diag(S);
if min(s)<=n*eps(max(s))
    error('polar_unitary:Singular','A must be numerically nonsingular.');
end
U=L*R';
end
