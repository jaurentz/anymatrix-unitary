function alpha = js_inverse_geronimus(diagJ,subJ,count)
%JS_INVERSE_GERONIMUS Recover real Verblunsky coefficients from Jacobi data.
%   ALPHA = JS_INVERSE_GERONIMUS(DIAGJ,SUBJ,K) uses the endpoint
%   Riccati form of the inverse Geronimus relations on [-2,2].
%   Requires at least floor(K/2)+1 diagonal and floor(K/2) subdiagonal
%   Jacobi coefficients. This is an independent check against the
%   explicit Jacobi-Szego Verblunsky formula.
validateattributes(count,{'numeric'},{'scalar','integer','nonnegative'});
alpha = zeros(count,1);
if count == 0, return; end
m = floor(count/2);
assert(numel(diagJ)>=m+1 && numel(subJ)>=m, ...
    'Insufficient Jacobi recurrence coefficients.');
assert(isreal(diagJ) && isreal(subJ) && all(subJ(:)>0), ...
    'Real Jacobi data with positive offdiagonal required.');
u = 2-diagJ(1);
v = 2+diagJ(1);
for j = 0:m
    if j >= 1
        assert(u>0 && v>0,'Riccati denominator is not positive.');
        u = 2-diagJ(j+1)-subJ(j)^2/u;
        v = 2+diagJ(j+1)-subJ(j)^2/v;
        if 2*j <= count
            alpha(2*j) = 1-(u+v)/2; % alpha_(2*j-1)
        end
    end
    assert(u>0 && v>0,'Jacobi data incompatible with support [-2,2].');
    if 2*j+1 <= count
        alpha(2*j+1) = (v-u)/(v+u); % alpha_(2*j)
    end
end
end
