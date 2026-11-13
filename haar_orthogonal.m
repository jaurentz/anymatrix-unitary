function [U,properties] = haar_orthogonal(n,seed,output)
%HAAR_ORTHOGONAL Haar-distributed random real orthogonal matrix.
%   U=HAAR_ORTHOGONAL(N,SEED) uses phase/sign-corrected Gaussian QR.
%   A specified integer SEED restores the caller's RNG state on exit.
%   Reference: F. Mezzadri, Notices AMS 54 (2007), 592-604.
%   https://arxiv.org/abs/math-ph/0609050
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [14] https://doi.org/10.48550/arXiv.math-ph/0609050 (arXiv preprint DOI)
properties={'square','unitary','normal','orthogonal','real','random','scalable'};
% Optional analytic spectral-reference mode.
if nargin>=3
    if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
        U=unitary_spectral_reference('haar_orthogonal',n,seed);
        return
    end
    error('haar_orthogonal:Mode','Only terminal spectral mode is supported for this generator.');
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin>=2 && ~isempty(seed)
    validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','finite'});
    old=rng; restore=onCleanup(@()rng(old)); %#ok<NASGU>
    rng(seed,'twister');
end
[Q,R]=qr(randn(n));
d=diag(R); s=ones(n,1); s(d<0)=-1;
U=Q*diag(s);
end
