function [U, properties] = cycle_propagator(n, t, output)
%CYCLE_PROPAGATOR Continuous-time quantum propagator on a periodic chain.
%   U = CYCLE_PROPAGATOR(N,T) returns exp(-1i*T*H), where
%   H = P + P', and P is the N-by-N cyclic shift.
%   T defaults to 1. N=1,2 use this same algebraic definition of H,
%   which counts the periodic bonds with their natural multiplicities.
%
%   The exact eigenvalues are exp(-2i*T*cos(2*pi*k/N)), k=0:N-1.
%   We build U from the inverse FFT of these eigenvalues, without EXPM.
%   Its limiting phase law is the pushforward of the uniform circle law
%   under theta -> -2*T*cos(theta) modulo 2*pi; for larger |T| the
%   phase map wraps around the unit circle.
%
%   Reference:
%   N. W. Ashcroft and N. D. Mermin, Solid State Physics, 1976
%   (periodic tight-binding dispersion relation).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [2] https://books.google.com/books?id=1C9HAQAAIAAJ
properties = {'unitary', 'square', 'circulant', 'scalable', ...
              'complex', 'eigenvalues', 'parameter-dependent'};
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
if nargin < 2
    t = 1;
end
validateattributes(t, {'numeric'}, {'scalar', 'real', 'finite'});
k = (0:n-1)';
lambda = 2*cos(2*pi*k/n);
z = exp(-1i*t*lambda);
if nargin<3,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_fft_apply(X,flag,'circulant',z,[]);
    return
end
c = ifft(z);
idx = mod((0:n-1)' - (0:n-1), n) + 1;
U = c(idx);
end
