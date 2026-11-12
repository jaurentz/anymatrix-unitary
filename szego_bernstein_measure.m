function mu=szego_bernstein_measure(a)
%SZEGO_BERNSTEIN_MEASURE Explicit Bernstein-Szego circle weight.
%   MU=SZEGO_BERNSTEIN_MEASURE(A) gives the normalized Poisson weight
%   (1-|A|^2)/(2*pi*|exp(i*theta)-conj(A)|^2) on the unit circle.
%   The Verblunsky sequence is A,0,0,... with |A|<1.
if nargin<1,a=0.4;end
mu=js_bernstein_measure(a);
end
