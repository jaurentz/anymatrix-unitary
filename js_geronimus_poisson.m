function rho=js_geronimus_poisson(theta,a,r)
%JS_GERONIMUS_POISSON Exact Poisson-smoothed infinite Geronimus measure.
%   RHO=JS_GERONIMUS_POISSON(THETA,A,R) evaluates
%       (P_R * mu_a)(theta) = Re F(R*exp(i*theta))/(2*pi)
%   for the infinite constant-Verblunsky CMV e1 measure. This includes
%   the atom when present. For a=0 the result is uniform.
%   This is NOT smoothing the equal-weight finite trace measure.
%
%   The Schur fixed point solves
%       conj(a)*z*f(z)^2 + (1-z)*f(z) - a = 0,
%   and the analytic branch at z=0 is
%       f(z) = 2*a/(1-z+sqrt((1-z)^2+4*abs(a)^2*z)).
%   F(z)=(1+z*f(z))/(1-z*f(z)).
%
%   M. J. Colbrook, Commun. Math. Phys. 384 (2021), Section 6.3:
%   https://doi.org/10.1007/s00220-021-04072-4
if nargin<3,error('js_geronimus_poisson:Arguments','Supply theta, a, and r.');end
validateattributes(theta,{'numeric'},{'real','finite'});
validateattributes(a,{'numeric'},{'scalar','finite'});
validateattributes(r,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
if abs(a)>=1
    error('js_geronimus_poisson:Alpha','A must lie in the open disk.');
end
z=r*exp(1i*theta);
if a==0
    rho=ones(size(theta))/(2*pi);
else
    discriminant=(1-z).^2+4*abs(a)^2*z;
    f=(2*a)./(1-z+sqrt(discriminant));
    g=z.*f;
    rho=max(0,real((1+g)./(1-g))/(2*pi));
end
end
