function w = js_jacobi_circle_density(theta,a,b)
%JS_JACOBI_CIRCLE_DENSITY Normalized Szego lift density with respect to dtheta.
%   w(theta) = (1-cos(theta))^(a+1/2)*(1+cos(theta))^(b+1/2)
%              / (2^(a+b+2)*B(a+1,b+1)); a,b>-1.
%   Integrates to 1 over [-pi,pi]. Endpoint values may be Inf if the
%   corresponding exponent is negative (integrable singularities).
logB = gammaln(a+1) + gammaln(b+1) - gammaln(a+b+2);
t = cos(theta);
w = (1-t).^(a+0.5) .* (1+t).^(b+0.5) .* exp(-(a+b+2)*log(2)-logB);
end
