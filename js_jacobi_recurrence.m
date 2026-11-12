function [diagJ,subJ] = js_jacobi_recurrence(a,b,m)
%JS_JACOBI_RECURRENCE Jacobi recurrence for x=2*t, x in [-2,2].
%   [DIAGJ,SUBJ] returns m diagonal and m-1 positive subdiagonal
%   coefficients in the *orthonormal* Jacobi three-term recurrence.
%   The first off-diagonal is evaluated separately to avoid a removable
%   0/0 when a+b=-1. Independent of the Verblunsky formulas.
validateattributes(a,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(b,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(m,{'numeric'},{'scalar','integer','nonnegative'});
diagJ = zeros(m,1);
subJ = zeros(max(0,m-1),1);
if m == 0, return; end
diagJ(1) = 2*(b-a)/(a+b+2);
for j = 1:m-1
    diagJ(j+1) = 2*(b^2-a^2)/((2*j+a+b)*(2*j+a+b+2));
end
if m >= 2
    subJ(1) = sqrt(16*(a+1)*(b+1)/((a+b+2)^2*(a+b+3)));
end
for j = 1:m-2
    subJ(j+1) = sqrt(16*(j+1)*(j+a+b+1)*(j+a+1)*(j+b+1)/ ...
        ((2*j+a+b+1)*(2*j+a+b+2)^2*(2*j+a+b+3)));
end
end
