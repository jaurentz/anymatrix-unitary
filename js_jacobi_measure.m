function mu = js_jacobi_measure(a,b)
%JS_JACOBI_MEASURE Normalized Jacobi measure and its Szego circle lift.
%   MU = JS_JACOBI_MEASURE(A,B), A,B > -1, describes the Jacobi
%   probability measure (1-t)^A(1+t)^B dt on [-1,1], its image on
%   [-2,2] by x=2*t, and its conjugation-symmetric Szego lift on T.
%   MU.CIRCLE_DENSITY(THETA) is the density with respect to dtheta,
%   not dtheta/(2*pi); similarly MU.REAL_DENSITY(X) is w.r.t. dx.
if nargin < 1 || isempty(a), a = 0; end
if nargin < 2 || isempty(b), b = 0; end
validateattributes(a,{'numeric'},{'real','finite','scalar','>',-1});
validateattributes(b,{'numeric'},{'real','finite','scalar','>',-1});
mu = struct();
mu.family = 'jacobi_szego';
mu.measure_type = 'absolutely_continuous';
mu.a = a;
mu.b = b;
mu.real_support = [-2 2];
mu.unit_circle_support = 'full circle (closure), possible endpoint zeros/singularities';
mu.map = 'x = z + z^(-1) = 2*cos(theta)';
mu.real_density = @(x) js_jacobi_real_density(x,a,b);
mu.circle_density = @(theta) js_jacobi_circle_density(theta,a,b);
mu.log_beta = gammaln(a+1) + gammaln(b+1) - gammaln(a+b+2);
end
