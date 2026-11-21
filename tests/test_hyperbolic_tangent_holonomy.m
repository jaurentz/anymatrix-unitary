function tests=test_hyperbolic_tangent_holonomy
%TEST_HYPERBOLIC_TANGENT_HOLONOMY Intrinsic H^n parallel transport.
tests=functiontests(localfunctions);
end

function test_exact_triangle_area_and_Poincare_disk(tc)
for n=[2 3 4 7]
    pts=triangle(n);
    U=hyperbolic_tangent_holonomy(n,pts);
    S=hyperbolic_tangent_holonomy(n,pts,'spectral');
    E=area_independent(pts);
    verifyLessThan(tc,abs(S.triangle_area-E),2e-11);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.dimension,n);
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),2e-11*n);
    verifyEqual(tc,det(U),1,'AbsTol',2e-11*n);
    lambda=S.finite.eigenvalues();
    verifyEqual(tc,numel(lambda),n);
    vals=eig(U);
    verifyLessThan(tc,max(arrayfun(@(x)min(abs(vals-x)),lambda)), ...
       4e-10);
    for k=-3:5
        verifyLessThan(tc,abs(trace(U^abs(k))/n ...
            -S.finite.moments(abs(k))),2e-10);
    end
    verifyLessThan(tc,norm(U-eye(n),'fro'), ...
       2*sqrt(2)*abs(sin(E/2))+1e-10);
    if n==2
        z=(pts(2,:)+1i*pts(3,:))./(1+pts(1,:));
        baseline=hyperbolic_triangle_holonomy(2,z);
        verifyLessThan(tc,norm(U-baseline,'fro'),2e-11);
    end
end
end

function test_nonplanar_H3_and_two_noncommuting_holonomies(tc)
n=3;home=[1;zeros(n,1)];
a=point([.59;.0;.0]);b=point([0;.57;0]);
c=point([0;0;.47]);
L1=[home,a,b];L2=[home,b,c];
U1=hyperbolic_tangent_holonomy(n,L1);
U2=hyperbolic_tangent_holonomy(n,L2);
verifyGreaterThan(tc,norm(U1*U2-U2*U1,'fro'),5e-3);
quad=[home,a,b,c];
D=hyperbolic_tangent_holonomy(3,quad,'model');
verifyLessThan(tc,norm(D.U'*D.U-eye(3),'fro'),4e-11);
verifyLessThan(tc,D.orientation_determinant-1,3e-11);
verifyGreaterThan(tc,norm(D.U-eye(3),'fro'),.1);
J=diag([-1,ones(1,n)]);
verifyLessThan(tc,norm(D.lorentz'*J*D.lorentz-J,'fro'), ...
    1e-10);
verifyLessThan(tc,norm(D.lorentz*home-home),1e-10);
verifyEqual(tc,det(D.lorentz),1,'AbsTol',1e-10);
for i=1:size(quad,2)
    B=D.leg_boosts(:,:,i);
    verifyLessThan(tc,norm(B'*J*B-J,'fro'),2e-10);
    verifyEqual(tc,det(B),1,'AbsTol',2e-10);
    verifyLessThan(tc,norm(B*quad(:,i)- ...
        quad(:,mod(i,size(quad,2))+1)),3e-10);
end
S=hyperbolic_tangent_holonomy(3,quad,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
end

function test_reverse_path_and_tangent_gauge_covariance(tc)
for n=[2 3 5]
    pts=triangle(n);
    loop=[pts(:,1),pts(:,3),pts(:,2)];
    U=hyperbolic_tangent_holonomy(n,pts);
    V=hyperbolic_tangent_holonomy(n,loop);
    verifyLessThan(tc,norm(U*V-eye(n),'fro'),3e-10);
    D=hyperbolic_tangent_holonomy(n,pts,'model');
    steps=D.leg_tangent_maps;
    accumulated=eye(n);
    for i=1:size(steps,3)
        accumulated=steps(:,:,i)*accumulated;
    end
    verifyLessThan(tc,norm(accumulated-U,'fro'),3e-10);
    % Replace EACH of the canonical orthonormal frames by any
    % independent SO(n) gauge rotation. The loop changes by
    % similarity in the basepoint gauge, not arbitrary entries.
    gauges=cell(3,1);
    for i=1:3
        A=zeros(n);
        for k=1:n-1,A(k,k+1)=.16*sin(k+.33*i);end
        gauges{i}=expm(A-A');
    end
    change=eye(n);
    for i=1:3
        next=mod(i,3)+1;
        R=gauges{next}'*steps(:,:,i)*gauges{i};
        change=R*change;
    end
    verifyLessThan(tc,norm(change-gauges{1}'*U*gauges{1},'fro'), ...
        3e-10);
end
end

function test_geodesic_velocity_and_independent_connection_ODE(tc)
for n=[2 3 5]
    pts=triangle(n);J=diag([-1,ones(1,n)]);
    for i=1:3
        x=pts(:,i);y=pts(:,mod(i,3)+1);
        [B,R,Ex,P]=um_hyperbolic_leg(x,y,J);
        Ey=um_hyperboloid_frame(y);
        verifyLessThan(tc,norm(x'*J*Ex),2e-12);
        verifyLessThan(tc,norm(Ex'*J*Ex-eye(n),'fro'),2e-12);
        verifyLessThan(tc,norm(y'*J*P),2e-12);
        verifyLessThan(tc,norm(P'*J*P-eye(n),'fro'),2e-12);
        verifyLessThan(tc,norm(R'*R-eye(n),'fro'),2e-12);
        verifyLessThan(tc,norm(B*x-y),2e-12);
        O=integrate_connection_RK4(x,y,Ex,300);
        verifyLessThan(tc,norm(P-O,'fro'),3e-8);
        c=-x'*J*y;dist=acosh(c);
        tangent=(y-c*x)/sinh(dist);
        out=tangent+(x+y)*(y'*J*tangent)/(1+c);
        tangentEnd=(c*y-x)/sinh(dist);
        verifyLessThan(tc,norm(out-tangentEnd),3e-11);
    end
end
end

function test_complex_polygon_higher_dimension_and_rejection(tc)
for n=[3 4 8]
    U=hyperbolic_tangent_holonomy(n);
    S=hyperbolic_tangent_holonomy(n,[],'spectral');
    verifyEqual(tc,size(U),[n n]);
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),5e-10);
    verifyEqual(tc,det(U),1,'AbsTol',5e-10);
    verifyEqual(tc,S.finite.kind,'unknown');
end
[A,prop]=hyperbolic_tangent_holonomy();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(prop,'orthogonal')));
verifyTrue(tc,any(strcmp(prop,'unitary')));
verifyError(tc,@()hyperbolic_tangent_holonomy(1), ...
    'hyperbolic_tangent_holonomy:Order');
verifyError(tc,@()hyperbolic_tangent_holonomy(3,rand(4,2)), ...
    'hyperbolic_tangent_holonomy:Points');
P=triangle(3);
P(1,2)=P(1,2)+.03;
verifyError(tc,@()hyperbolic_tangent_holonomy(3,P), ...
    'hyperbolic_tangent_holonomy:Points');
P=triangle(3);P(:,2)=P(:,1);
verifyError(tc,@()hyperbolic_tangent_holonomy(3,P), ...
    'hyperbolic_tangent_holonomy:Leg');
verifyError(tc,@()hyperbolic_tangent_holonomy(3,[], ...
    'nonsense'),'hyperbolic_tangent_holonomy:Mode');
end

function pts=triangle(n)
p=zeros(n,1);q=zeros(n,1);
p(1)=.51;p(2)=.12;
q(1)=.10;q(2)=.47;
if n>=3,p(3)=.18;q(3)=-.11;end
pts=[[1;zeros(n,1)],point(p),point(q)];
end

function x=point(v)
rho=norm(v);
x=[cosh(rho);sinh(rho)*v/rho];
end

function area=area_independent(pts)
% Independent hyperbolic L'Huilier formula, not corner-angle code.
J=diag([-1,ones(1,size(pts,1)-1)]);
a=acosh(-pts(:,2)'*J*pts(:,3));
b=acosh(-pts(:,3)'*J*pts(:,1));
c=acosh(-pts(:,1)'*J*pts(:,2));
s=(a+b+c)/2;
area=4*atan(sqrt(tanh(s/2)*tanh((s-a)/2) ...
    *tanh((s-b)/2)*tanh((s-c)/2)));
end

function V=integrate_connection_RK4(x,y,V,steps)
% gamma(t) joins x->y on hyperboloid. Parallel equation in R^(n+1)
% V'=<gamma',V>_Lorentz*gamma; integrate with separate RK4.
J=diag([-1,ones(1,numel(x)-1)]);
c=-x'*J*y;dist=acosh(c);
rhs=@(t,W)rhs_eval(t,W,x,y,dist,J);
h=1/steps;
for j=1:steps
    t=(j-1)*h;
    k1=rhs(t,V);
    k2=rhs(t+h/2,V+h*k1/2);
    k3=rhs(t+h/2,V+h*k2/2);
    k4=rhs(t+h,V+h*k3);
    V=V+(h/6)*(k1+2*k2+2*k3+k4);
end
end

function dV=rhs_eval(t,V,x,y,a,J)
g=sinh((1-t)*a)/sinh(a)*x+sinh(t*a)/sinh(a)*y;
dg=-a*cosh((1-t)*a)/sinh(a)*x ...
    +a*cosh(t*a)/sinh(a)*y;
dV=g*(dg'*J*V);
end
