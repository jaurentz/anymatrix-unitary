function [U,properties] = circular_beta(n,dysonBeta,seed,output)
%CIRCULAR_BETA Sparse random CMV model for Dyson's circular beta ensemble.
%   U=CIRCULAR_BETA(N,BETA,SEED), BETA>0, generates the Killip-Nenciu
%   Verblunsky coefficient ensemble and closes it with a uniform phase.
%   The eigenphase joint density is proportional to
%   prod_(j<k) |exp(1i*theta_j)-exp(1i*theta_k)|^BETA.
%   Dyson BETA=2 reproduces the CUE eigenvalue distribution (but the
%   matrix is sparse CMV, not Haar-distributed as an N-by-N matrix).
%   BETA=1,4 correspond to circular beta eigenvalue laws with their
%   own eigenvalue counting conventions, not Haar O(N) or USp(N).
%   Reference: R. Killip and I. Nenciu, IMRN 2004(50), 2665-2701.
%   https://doi.org/10.1155/S1073792804141597
properties={'square','unitary','normal','sparse','banded','random', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
if nargin<2 || isempty(dysonBeta),dysonBeta=2;end
if nargin<3,seed=[];end
[alpha,betaTerminal]=js_circular_beta_coeffs(n,dysonBeta,seed);
if nargin<4,output='matrix';end
U=js_core_product(alpha,betaTerminal,'cmv',[],output);
end
