function mu = szego_geronimus_measure(a)
%SZEGO_GERONIMUS_MEASURE OPUC measure with constant Verblunsky parameter.
%   MU=SZEGO_GERONIMUS_MEASURE(A), |A|<1, has alpha_k=A for all k>=0.
%   Its spectral measure can include an absolutely continuous arc and
%   a point mass. Analytic continuous arc density, atom and exact
%   Poisson smoothing are exposed as distinguished-vector metadata.
%   No universal finite trace-density statement is implied.
%   Reference: Simon, OPUC (2005), Example 1.6.12.
%   Primary reference DOI / stable reference:
%   https://doi.org/10.1090/coll054.1
if nargin<1 || isempty(a),a=0.3+0.2i;end
validateattributes(a,{'numeric'},{'scalar','finite'});
if abs(a)>=1,error('szego_geronimus_measure:Alpha','A must lie in the open disk.');end
% Analytic decomposition of the INFINITE-dimensional e1 measure.
% Colbrook, Commun. Math. Phys. 384 (2021), Section 6.3:
% https://doi.org/10.1007/s00220-021-04072-4
% Atom at b=2*arg(1+conj(a)) occurs iff real(a)+abs(a)^2>0.
% Its mass is max(0,2*(real(a)+abs(a)^2)/abs(1+a)^2).
% Neither density nor atom is the finite CMV equal-weight trace measure.
edge=2*asin(abs(a));
b=2*angle(1+conj(a));
mass=max(0,2*(real(a)+abs(a)^2)/abs(1+a)^2);
if mass>0
    atoms=struct('angles',b,'weights',mass);
    type='mixed';
else
    atoms=struct('angles',[],'weights',[]);
    type='absolutely_continuous';
end
if a==0,support=[-pi pi];else,support=[-pi -edge;edge pi];end
mu=struct('family','geronimus','measure_type',type, ...
    'alpha',a,'verblunsky_class','real_or_complex', ...
    'unit_circle_support','parameter-dependent essential arc', ...
    'ac_support',support,'point_mass',atoms, ...
    'ac_mass',1-mass, ...
    'circle_density',@(theta)js_geronimus_density(theta,a), ...
    'poisson_density',@(theta,r)js_geronimus_poisson(theta,a,r));
end
