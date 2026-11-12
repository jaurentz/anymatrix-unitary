function mu = js_bernstein_measure(a)
%JS_BERNSTEIN_MEASURE A one-parameter Bernstein-Szego/Poisson measure.
%   For a with |a|<1 the normalized circle measure (w.r.t. dtheta) is
%       (1-|a|^2)/(2*pi*|exp(i*theta)-conj(a)|^2).
%   Its Verblunsky coefficients are alpha_0=a, alpha_j=0 for j>=1.
%   This convention gives integral z dmu(z)=conj(a).
if nargin<1 || isempty(a),a=0.4;end
validateattributes(a,{'numeric'},{'scalar','finite'});
if abs(a)>=1,error('js_bernstein_measure:Alpha','The Poisson parameter must have modulus < 1.');end
mu=struct();
mu.family='bernstein_szego';
mu.measure_type='absolutely_continuous';
mu.alpha0=a;
mu.circle_density=@(theta) (1-abs(a)^2)./(2*pi*abs(exp(1i*theta)-conj(a)).^2);
mu.unit_circle_support='full circle';
end
