function [B,R,Ex,P]=um_hyperbolic_leg(x,y,J)
%UM_HYPERBOLIC_LEG Intrinsic Levi-Civita PT and its Lorentz boost.
%   P=PT_x->y(Ex) has tangent Gram I at y.
%   B=[y,P]*J*[x,Ex]'*J belongs to SO+(1,n), sends x to y.
n=numel(x)-1;c=-x'*J*y;
Ex=um_hyperboloid_frame(x);
Ey=um_hyperboloid_frame(y);
P=Ex+(x+y)*(y'*J*Ex)/(1+c);
Qx=[x,Ex];Qy=[y,P];
B=Qy*J*Qx'*J;
R=Ey'*J*P;
if norm(R'*R-eye(n),'fro')>5e-8*n ...
        || norm(B'*J*B-J,'fro')>5e-8*n
    error('um_hyperbolic_leg:Numerics', ...
        'Tangent transport failed metric preservation.');
end
end
