function [U,alpha] = js_measure_to_matrix(mu,n,representation,beta,twists,output)
%JS_MEASURE_TO_MATRIX Generic measure -> Verblunsky -> finite core matrix.
%   [U,ALPHA]=JS_MEASURE_TO_MATRIX(MU,N,REP,BETA,TWISTS)
%   takes an OPUC probability measure descriptor, generates the first
%   N-1 interior Verblunsky coefficients, and closes at unit phase BETA.
%   Rep and twists are passed to JS_CORE_PRODUCT.
if nargin<3 || isempty(representation), representation='cmv'; end
if nargin<4 || isempty(beta), beta=1; end
if nargin<5, twists=[]; end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
alpha=js_verblunsky(mu,n-1);
if nargin<6,output='matrix';end
U=js_core_product(alpha,beta,representation,twists,output);
end
