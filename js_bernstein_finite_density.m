function w=js_bernstein_finite_density(theta,alpha)
%JS_BERNSTEIN_FINITE_DENSITY Rational Bernstein-Szego density wrt dtheta.
z=exp(1i*theta);
p=ones(size(z));pstar=p;
for j=1:numel(alpha)
    old=p;
    p=z.*p-conj(alpha(j)).*pstar;
    pstar=pstar-alpha(j).*z.*old;
end
c=prod(1-abs(alpha(:)).^2);
w=c./(2*pi*abs(pstar).^2);
end
