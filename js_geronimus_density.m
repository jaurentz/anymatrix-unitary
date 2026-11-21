function rho=js_geronimus_density(theta,a)
%JS_GERONIMUS_DENSITY Absolutely continuous infinite OPUC e1 density.
%   RHO=JS_GERONIMUS_DENSITY(THETA,A) returns dmu_ac/dtheta for the
%   infinite CMV operator with constant Verblunsky alpha_j=A, |A|<1.
%   THETA is in radians (periodic). The full probability measure may
%   ALSO contain a point mass: see SZEGO_GERONIMUS_MEASURE.
%
%   Let edge=2*asin(abs(A)) and b=2*arg(1+conj(A)). On the
%   essential arc |theta|>edge, with theta wrapped to [-pi,pi],
%
%   rho(theta) = sqrt(sin(theta/2)^2-|A|^2)
%       / (2*pi*abs(1+A)*abs(sin((theta-b)/2))).
%
%   This is an INFINITE-dimensional distinguished-vector measure,
%   NOT the finite trace spectral density of GERONIMUS_CMV.
%   a=0 is exactly uniform with rho=1/(2*pi).
%
%   M. J. Colbrook (2021), Section 6.3, formula for Geronimus measure:
%   https://doi.org/10.1007/s00220-021-04072-4
%   B. Simon, OPUC (2005): https://doi.org/10.1090/coll054.1
if nargin<2,error('js_geronimus_density:Arguments','Supply theta and a.');end
validateattributes(theta,{'numeric'},{'real','finite'});
validateattributes(a,{'numeric'},{'scalar','finite'});
if abs(a)>=1
    error('js_geronimus_density:Alpha','A must lie in the open disk.');
end
if a==0
    rho=ones(size(theta))/(2*pi);
    return
end
phase=mod(theta+pi,2*pi)-pi;
edge=2*asin(abs(a));
b=2*angle(1+conj(a));
rho=zeros(size(theta));
inside=abs(phase)>edge;
if any(inside(:))
    p=phase(inside);
    top=sqrt(max(0,sin(p/2).^2-abs(a)^2));
    rho(inside)=top./(2*pi*abs(1+a)*abs(sin((p-b)/2)));
end
end
