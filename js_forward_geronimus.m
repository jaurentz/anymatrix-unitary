function [diagJ,subJ] = js_forward_geronimus(alpha,m)
%JS_FORWARD_GERONIMUS Jacobi recurrence from real Verblunsky coefficients.
%   Uses the direct Geronimus relations on [-2,2], alpha_-1=-1.
%   Requires alpha_0,...,alpha_(2*m-1). Output has m diagonal and
%   m positive subdiagonal values, including the coupling beyond an
%   m-by-m principal Jacobi section.
%   Primary reference DOI / stable reference:
%   https://doi.org/10.1090/S0273-0979-05-01075-X
validateattributes(m,{'numeric'},{'scalar','integer','nonnegative'});
assert(isreal(alpha) && all(abs(alpha(:))<1), ...
    'Forward Geronimus requires real interior Verblunsky coefficients.');
assert(numel(alpha) >= 2*m,'Need at least 2*m coefficients.');
diagJ = zeros(m,1);
subJ = zeros(m,1);
for j = 0:m-1
    if j == 0
        am1 = -1; am2 = 0;
    else
        am1 = alpha(2*j);   % alpha_(2*j-1)
        am2 = alpha(2*j-1); % alpha_(2*j-2)
    end
    aeven = alpha(2*j+1);
    aodd = alpha(2*j+2);
    diagJ(j+1) = (1-am1)*aeven - (1+am1)*am2;
    subJ(j+1) = sqrt((1-am1)*(1-aeven^2)*(1+aodd));
end
end
