function [U,properties]=hyperbolic_tangent_holonomy(n,points,mode)
%HYPERBOLIC_TANGENT_HOLONOMY Levi-Civita polygon holonomy in H^n.
%   U=HYPERBOLIC_TANGENT_HOLONOMY(N,POINTS,MODE) is the genuine
%   SO(N) parallel-transport holonomy of a closed geodesic polygon on
%       H^n={x in R^(n+1): -x(1)^2+sum_j x(j+1)^2=-1, x(1)>0}.
%   POINTS has n+1 rows and at least three columns of Lorentz-unit
%   future-pointing vertices; the closing last->first edge is implicit.
%   Points need not be coplanar. N>=2; for N>=3 the default four-point
%   loop uses 3 spatial coordinates (genuinely noncoplanar in H^3).
%   N=2 default uses a nondegenerate geodesic triangle.
%
%   In Lorentz metric J=diag(-1,I), each geodesic leg x->y has exact
%   parallel transport of a tangent vector V (x'*J*V=0):
%       PT(V)=V+(x+y)*(y'*J*V)/(1-x'*J*y).
%   Build canonical, Lorentz-orthonormal tangent frames
%       E(x)=[x_spatial'; I+x_spatial*x_spatial'/(1+x(1))]
%   and the uniquely determined proper-orthochronous Lorentz boost
%       B_xy=[y,PT(E(x))]*J*[x,E(x)]'*J.
%   It maps x to y and all tangent vectors by Levi-Civita PT.
%   For a polygon with base point x0, the Lorentz product B_loop fixes
%   x0, and U=E(x0)'*J*B_loop*E(x0) is in SO(N). The ambient B_loop
%   is in SO^+(1,N) and is NOT Euclidean orthogonal/unitary.
%
%   MODE='matrix' (default) returns the n-by-n real orthogonal
%   tangent-fiber holonomy; MODE='model' exposes boosts, frames,
%   Lorentz total product, vertex/leg transports and returned U.
%   MODE='spectral': triangles have EXACT finite eigenphase law
%   exp(+-i hyperbolic defect area), with n-2 unit eigenvalues;
%   arbitrary multi-plane polygons have UNKNOWN closed analytic finite
%   phase law (do not use eig(U) as a purported exact phase reference).
%
%   This is distinct from hyperbolic_triangle_holonomy, whose H^2
%   angle is merely identity-padded to higher matrix order.
%
%   [142] Nagano et al., A Wrapped Normal Distribution on Hyperbolic Space,
%   ICML (2019), source for hyperboloid geodesics / tangent PT:
%   https://proceedings.mlr.press/v97/nagano19a.html
%   [113] Hyperbolic geodesic holonomy baseline (Poincare disk).
%   https://doi.org/10.1515/9780691219899
properties={'square','unitary','orthogonal','real','normal', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
      &&n==floor(n)&&n>=2&&n<=40)
    error('hyperbolic_tangent_holonomy:Order','N must be 2..40.');
end
if nargin<2||isempty(points)
    P=zeros(n+1,4);
    if n==2,P=P(:,1:3);end
    P(1,1)=1;
    if n==2
        w=[.51,.1;.12,.47];
    else
        w=[.52,-.16,.25;.1,.43,-.11;.14,.18,.39];
    end
    for k=2:size(P,2)
        vec=zeros(n,1);
        vec(1:size(w,1))=w(:,k-1);
        rho=norm(vec);
        P(:,k)=[cosh(rho);sinh(rho)*vec/rho];
    end
    points=P;
end
if nargin<3||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('hyperbolic_tangent_holonomy:Mode','MODE must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','spectral'}))
    error('hyperbolic_tangent_holonomy:Mode','Unknown mode.');
end
if ~(isnumeric(points)&&isreal(points)&&ismatrix(points) ...
     &&size(points,1)==n+1&&size(points,2)>=3 ...
     &&size(points,2)<=24&&all(isfinite(points(:))) ...
     &&all(points(1,:)>0)&&max(abs(points(:)))<1e4)
    error('hyperbolic_tangent_holonomy:Points', ...
        'POINTS must contain 3..24 finite future hyperboloid vertices.');
end
m=size(points,2);
J=diag([-1,ones(1,n)]);
for j=1:m
    x=points(:,j);
    if abs(x'*J*x+1)>3e-9*max(1,x(1)^2)
        error('hyperbolic_tangent_holonomy:Points', ...
            'Every vertex must have Lorentz norm -1.');
    end
    y=points(:,mod(j,m)+1);
    c=-x'*J*y;
    if c<=1+1e-11
        error('hyperbolic_tangent_holonomy:Leg', ...
            'Repeated, numerically degenerate or invalid leg.');
    end
end
if strcmp(mode,'spectral')
    empty=struct('kind','unknown','density',[],'support',[], ...
        'atoms',[],'moments',[]);
    finite=struct('kind','unknown','eigenvalues',[], ...
        'moments',[],'atoms',[]);
    notes=['Generic polygon spanning several totally geodesic ', ...
        'planes in H^n can yield non-Abelian SO(n) holonomy; ', ...
        'no generic analytic trace eigenphase law is asserted. ', ...
        'Lorentz boosts themselves are not unitary matrices.'];
    if m==3
        area=um_hyperbolic_triangle_area(points,J);
        finite.kind='exact';
        finite.eigenvalues=@()[exp(1i*area);exp(-1i*area); ...
            ones(n-2,1)];
        finite.moments=@(k)(n-2+2*cos(k*area))/n;
        if n==2
            finite.atoms=struct('angles',[area;-area], ...
                'weights',[1;1]/n);
        else
            finite.atoms=struct('angles',[area;-area;0], ...
                'weights',[1;1;n-2]/n);
        end
        notes=['Every nondegenerate hyperboloid geodesic triangle lies ', ...
            'in a totally geodesic H^2. Exact defect area yields the ', ...
            'SO(2) rotation eigenpair with n-2 fixed tangent directions.'];
    end
    U=struct('family','hyperbolic_tangent_holonomy', ...
        'dimension',n, ...
        'normalization','probability (density with respect to dtheta)', ...
        'finite',finite,'limit',empty,'ensemble',empty,'underlying',empty, ...
        'reference','Nagano et al., ICML (2019), hyperboloid parallel transport', ...
        'notes',notes);
    if m==3,U.triangle_area=area;end
    return
end
B=eye(n+1);
if strcmp(mode,'model')
    legs=zeros(n+1,n+1,m);
    local=zeros(n,n,m);
    frames=zeros(n+1,n,m);
end
for j=1:m
    x=points(:,j);
    y=points(:,mod(j,m)+1);
    [step,localMap,Ex]=um_hyperbolic_leg(x,y,J);
    B=step*B;
    if strcmp(mode,'model')
        legs(:,:,j)=step;
        local(:,:,j)=localMap;
        frames(:,:,j)=Ex;
    end
end
E0=um_hyperboloid_frame(points(:,1));
U=real(E0'*J*B*E0);
if strcmp(mode,'model')
    U=struct('U',U,'lorentz',B,'lorentz_metric',J, ...
        'points',points,'frames',frames,'leg_boosts',legs, ...
        'leg_tangent_maps',local, ...
        'basepoint_error',norm(B*points(:,1)-points(:,1)), ...
        'orthogonality_error',norm(U'*U-eye(n),'fro'), ...
        'orientation_determinant',det(U), ...
        'notes',['The real SO(n) matrix U acts in an orthonormal ', ...
        'tangent frame at the base vertex. The Lorentz boosts ', ...
        'belong to SO^+(1,n), not Euclidean O(n+1).']);
end
end
