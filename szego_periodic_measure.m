function mu = szego_periodic_measure(pattern)
%SZEGO_PERIODIC_MEASURE Periodic OPUC Verblunsky sequence.
%   PATTERN is a nonempty vector strictly inside the open disk.
%   The resulting infinite-periodic CMV measure can have spectral bands
%   and gaps; the finite realization requires a terminal unit phase.
%   Reference: Simon, OPUC Part 2 (2005), periodic coefficients.
%   Primary reference DOI / stable reference:
%   https://doi.org/10.1090/coll/054.2
if nargin<1 || isempty(pattern),pattern=[0.25; 0.2i];end
validateattributes(pattern,{'numeric'},{'vector','nonempty','finite'});
if any(abs(pattern(:))>=1)
    error('szego_periodic_measure:Alpha','Pattern entries must lie in the open disk.');
end
mu=struct('family','periodic','measure_type','periodic_spectral', ...
    'pattern',pattern(:),'period',numel(pattern), ...
    'verblunsky_class','real_or_complex','unit_circle_support','spectral bands');
end
