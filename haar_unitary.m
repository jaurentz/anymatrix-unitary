function [U, properties] = haar_unitary(n, seed)
%HAAR_UNITARY Haar-distributed random unitary matrix from complex QR.
%   U = HAAR_UNITARY(N) uses the current MATLAB random stream.
%   U = HAAR_UNITARY(N,SEED) gives a reproducible draw using the
%   Mersenne Twister and restores the caller's original RNG state.
%
%   Given complex standard Gaussian Z = Q*R, multiply each column of Q
%   by the unit-modulus diagonal phase R(j,j)/abs(R(j,j)).
%   The resulting U is Haar-distributed on U(N) in exact arithmetic.
%   Its eigenvalues have CUE statistics; the expected one-point
%   eigenphase distribution is uniform, not the phases of one sample.
%
%   Reference:
%   F. Mezzadri, How to generate random matrices from the classical
%   compact groups, Notices Amer. Math. Soc. 54 (2007), 592-604.
%   https://arxiv.org/abs/math-ph/0609050
properties = {'unitary', 'square', 'complex', 'scalable', 'random'};
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
if nargin >= 2 && ~isempty(seed)
    validateattributes(seed, {'numeric'}, ...
                       {'scalar', 'integer', 'nonnegative', 'finite'});
    previous_rng_state = rng;
    restore_rng_state = onCleanup(@() rng(previous_rng_state)); %#ok<NASGU>
    rng(seed,'twister');
end
Z = (randn(n)+1i*randn(n))/sqrt(2);
[Q,R] = qr(Z);
d = diag(R);
phase = ones(n,1);
nonzero = abs(d) > 0;
phase(nonzero) = d(nonzero)./abs(d(nonzero));
U = Q*diag(phase);
end
