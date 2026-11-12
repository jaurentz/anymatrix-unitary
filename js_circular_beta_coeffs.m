function [alpha,betaTerminal] = js_circular_beta_coeffs(n,dysonBeta,seed)
%JS_CIRCULAR_BETA_COEFFS Killip-Nenciu finite random CMV coefficients.
%   Interior alpha_(j-1), j=1,...,n-1, are rotationally invariant with
%   |alpha_(j-1)|^2 ~ Beta(1,dysonBeta*(n-j)/2). Terminal beta is uniform
%   on the unit circle. All distributions generated using base RAND only.
%   Optional SEED restores the caller's RNG state.
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(dysonBeta,{'numeric'},{'scalar','real','positive','finite'});
if nargin>=3 && ~isempty(seed)
    validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','finite'});
    old=rng; restore=onCleanup(@()rng(old)); %#ok<NASGU>
    rng(seed,'twister');
end
alpha=zeros(n-1,1);
for j=1:n-1
    shape=dysonBeta*(n-j)/2;
    % Stable inverse CDF for a Beta(1,shape) radial square.
    radiusSquared=-expm1(log(rand)/shape);
    radiusSquared=min(1-8*eps,max(0,radiusSquared));
    alpha(j)=sqrt(radiusSquared)*exp(2*pi*1i*rand);
end
betaTerminal=exp(2*pi*1i*rand);
end
