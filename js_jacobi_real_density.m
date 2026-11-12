function w = js_jacobi_real_density(x,a,b)
%JS_JACOBI_REAL_DENSITY Jacobi probability density on [-2,2].
%   (2-x)^a (2+x)^b / (4^(a+b+1) * B(a+1,b+1)).
logB = gammaln(a+1) + gammaln(b+1) - gammaln(a+b+2);
w = zeros(size(x));
idx = abs(x) <= 2;
if any(idx(:))
    y = x(idx);
    w(idx) = (2-y).^a .* (2+y).^b .* exp(-2*(a+b+1)*log(2)-logB);
end
end
