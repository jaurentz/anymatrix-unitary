function alpha = js_verblunsky_jacobi(count,a,b)
%JS_VERBLUNSKY_JACOBI Explicit Verblunsky coefficients for a Jacobi Szego lift.
%   ALPHA = JS_VERBLUNSKY_JACOBI(K,A,B) returns alpha_0,...,alpha_(K-1).
%   alpha_(2*j)   = (b-a)/(2*j+a+b+2)
%   alpha_(2*j+1) = -(a+b+1)/(2*j+a+b+3).
%   K counts *interior* coefficients; a unit-modulus finite closure beta
%   is a separate independent parameter.
validateattributes(count,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(a,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(b,{'numeric'},{'scalar','real','finite','>',-1});
k = (0:count-1)';
alpha = zeros(count,1);
even = mod(k,2)==0;
alpha(even) = (b-a)./(k(even)+a+b+2);
alpha(~even) = -(a+b+1)./(k(~even)+a+b+2);
assert(all(abs(alpha)<1),'Jacobi parameters gave invalid Verblunsky coefficients');
end
