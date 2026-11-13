function mu = szego_geronimus_measure(a)
%SZEGO_GERONIMUS_MEASURE OPUC measure with constant Verblunsky parameter.
%   MU=SZEGO_GERONIMUS_MEASURE(A), |A|<1, has alpha_k=A for all k>=0.
%   Its spectral measure can include an absolutely continuous arc and
%   a point mass. No universal full-circle density is implied.
%   Reference: Simon, OPUC (2005), Example 1.6.12.
if nargin<1 || isempty(a),a=0.3+0.2i;end
validateattributes(a,{'numeric'},{'scalar','finite'});
if abs(a)>=1,error('szego_geronimus_measure:Alpha','A must lie in the open disk.');end
mu=struct('family','geronimus','measure_type','mixed_or_arc', ...
    'alpha',a,'verblunsky_class','real_or_complex', ...
    'unit_circle_support','parameter-dependent arc and possible atom');
end
