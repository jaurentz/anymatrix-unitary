function mu = szego_jacobi_measure(a,b)
%SZEGO_JACOBI_MEASURE Normalized Jacobi measure and symmetric Szego lift.
%   MU = SZEGO_JACOBI_MEASURE(A,B) describes the Jacobi probability
%   measure proportional to (1-t)^A*(1+t)^B on t in [-1,1] and its
%   conjugation-symmetric lift by x=2*t=z+z^(-1).
%   A,B > -1; by default A=B=0 (Legendre).
%
%   MU.REAL_DENSITY(x) is w.r.t. dx on [-2,2].
%   MU.CIRCLE_DENSITY(theta) is w.r.t. dtheta on [-pi,pi].
%   Use SZEGO_VERBLUNSKY(MU,K) to calculate coefficients.
%   Primary reference DOI / stable reference:
%   https://doi.org/10.1090/S0273-0979-05-01075-X
if nargin < 1, a=0; end
if nargin < 2, b=0; end
mu = js_jacobi_measure(a,b);
end
