function [angles,weights]=um_koopman_cyclic_measure(Q,H,c)
%UM_KOOPMAN_CYCLIC_MEASURE Numerical Koopman cyclic-vector phase weights.
%   c is a coefficient vector in the NONORTHOGONAL dictionary with
%   weighted metric G=H'*H. Transform to v=H*c, then compute weights
%   |V'*v|^2/(v'*v) in a unitary Schur basis of the whitened Q.
%   These are NUMERICAL FINITE spectral weights, not certified analytic
%   infinite-operator spectra. They generally differ from trace weights
%   1/n. Near-coincident Schur phases should be interpreted in clusters.
n=size(Q,1);
if ~isnumeric(c)||~isvector(c)||numel(c)~=n|| ...
   any(~isfinite(c(:)))||~any(c(:))
    error('um_koopman_cyclic_measure:Vector', ...
        'Coefficient vector must be nonzero and finite length N.');
end
v=H*c(:);
[W,T]=schur(Q,'complex');
angles=angle(diag(T));
weights=abs(W'*v).^2/real(v'*v);
end
