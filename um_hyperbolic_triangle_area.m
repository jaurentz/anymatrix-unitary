function area=um_hyperbolic_triangle_area(points,J)
%UM_HYPERBOLIC_TRIANGLE_AREA Positive hyperbolic angle defect.
%   Intrinsic H^2 area for any three hyperboloid vertices, even H^n.
%   Side a=dist(y,z), b=dist(z,x), c=dist(x,y).
x=points(:,1);y=points(:,2);z=points(:,3);
a=acosh(max(1,-y'*J*z));
b=acosh(max(1,-z'*J*x));
c=acosh(max(1,-x'*J*y));
ca=(cosh(b)*cosh(c)-cosh(a))/(sinh(b)*sinh(c));
cb=(cosh(c)*cosh(a)-cosh(b))/(sinh(c)*sinh(a));
cc=(cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b));
clip=@(v)min(1,max(-1,v));
area=pi-acos(clip(ca))-acos(clip(cb))-acos(clip(cc));
if ~isfinite(area)||area<=1e-12||area>=pi
    error('um_hyperbolic_triangle_area:Degenerate', ...
        'Require a noncollinear, well-resolved H^2 triangle.');
end
end
