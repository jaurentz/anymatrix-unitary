function E=um_hyperboloid_frame(x)
%UM_HYPERBOLOID_FRAME Canonical positive-spatial Lorentz tangent gauge.
%   x(1)>0 and -x(1)^2+norm(x(2:end))^2=-1.
%   J=diag(-1,I), E'*J*x=0 and E'*J*E=I.
v=x(2:end);n=numel(v);
E=[v';eye(n)+(v*v')/(1+x(1))];
end
