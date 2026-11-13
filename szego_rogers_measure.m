function mu = szego_rogers_measure(q)
%SZEGO_ROGERS_MEASURE Wrapped-normal Rogers-Szego probability measure.
%   0<Q<1. The wrapped Gaussian variance is -log(Q) and
%   alpha_k=(-1)^k*Q^((k+1)/2), k>=0 (real Verblunsky coefficients).
%   Reference: Simon, OPUC (2005), Example 1.6.5.
if nargin<1 || isempty(q),q=0.6;end
validateattributes(q,{'numeric'},{'scalar','real','finite','>',0,'<',1});
mu=struct('family','rogers_szego','measure_type','absolutely_continuous', ...
    'q',q,'verblunsky_class','real','unit_circle_support','full circle');
mu.circle_density=@(theta) js_rogers_density(theta,q);
end
