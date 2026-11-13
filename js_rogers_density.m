function w=js_rogers_density(theta,q)
%JS_ROGERS_DENSITY Wrapped-normal OPUC probability density wrt dtheta.
s=-log(q);
t=mod(theta+pi,2*pi)-pi;
M=ceil(1+sqrt(64*s)/(2*pi));
w=zeros(size(t));
for k=-M:M
    w=w+exp(-((t-2*pi*k).^2)/(2*s));
end
w=w/sqrt(2*pi*s);
end
