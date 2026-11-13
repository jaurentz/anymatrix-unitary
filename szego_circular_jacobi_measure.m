function mu = szego_circular_jacobi_measure(a,eta)
%SZEGO_CIRCULAR_JACOBI_MEASURE Positive generalized circular Jacobi weight.
%   MU=SZEGO_CIRCULAR_JACOBI_MEASURE(A,ETA), A>-1/2, ETA real, has
%   weight proportional to exp((pi-theta)*ETA)*|1-exp(i*theta)|^(2*A)
%   for theta in (0,2*pi). For ETA=0 the Verblunsky data are real.
%   See Simon, OPUC (2005); Bourgade, Nikeghbali, Rouault (2009).
if nargin<1 || isempty(a),a=0.3;end
if nargin<2 || isempty(eta),eta=0.5;end
validateattributes(a,{'numeric'},{'scalar','real','finite','>',-0.5});
validateattributes(eta,{'numeric'},{'scalar','real','finite'});
raw=@(t)exp((pi-t).*eta).*(4*sin(t/2).^2).^a;
normalization=integral(raw,0,pi,'RelTol',1e-11) + ...
    integral(raw,pi,2*pi,'RelTol',1e-11);
if ~isfinite(normalization) || normalization<=0
    error('szego_circular_jacobi_measure:Normalization', ...
        'Could not normalize the weight in double precision.');
end
mu=struct('family','circular_jacobi','measure_type','absolutely_continuous', ...
    'a',a,'eta',eta,'normalization',normalization, ...
    'verblunsky_class','real_or_complex', ...
    'unit_circle_support','full circle (closure)');
mu.circle_density=@(theta) raw(mod(theta,2*pi))./normalization;
end
