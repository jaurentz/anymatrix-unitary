function mu = szego_bernstein_finite_measure(a)
%SZEGO_BERNSTEIN_FINITE_MEASURE General finite Verblunsky perturbation.
%   For finite interior coefficient vector A (|A_j|<1), set later
%   coefficients to zero. This is a positive Bernstein-Szego weight
%   proportional to product(1-|A_j|^2)/|Phi_m^*(z)|^2.
%   A may be empty (uniform measure); use SZEGO_MATRIX for finite closure.
%   Reference: Simon, OPUC (2005); DLMF section 18.33.
if nargin<1,a=[0.3;0.2i];end
if ~isnumeric(a) || (~isempty(a) && ~isvector(a)) || any(~isfinite(a(:)))
    error('szego_bernstein_finite_measure:Alpha','A must be a finite numeric vector.');
end
if any(abs(a(:))>=1)
    error('szego_bernstein_finite_measure:Alpha','Every coefficient must lie in the open disk.');
end
mu=struct('family','bernstein_finite','measure_type','absolutely_continuous', ...
    'coefficients',a(:),'verblunsky_class','real_or_complex', ...
    'unit_circle_support','full circle');
mu.circle_density=@(theta) js_bernstein_finite_density(theta,a(:));
end
