function [H,properties]=geodesic_holonomy(n,points,output)
%GEODESIC_HOLONOMY Parallel transport around a spherical geodesic polygon.
%   H=GEODESIC_HOLONOMY(N,POINTS) uses real unit columns p_1,...,p_m
%   on S^(N-1), m>=3; implicit closing arc p_m -> p_1.
%   Each leg is the unique shortest geodesic between non-antipodal points.
%   H is the AMBIENT orthogonal product R(p_m,p_1)*...*R(p_1,p_2),
%   fixes p_1, and restricts to tangent holonomy at p_1.
%   For a triangle, all points lie in a 3-plane and the spectrum is
%   exp(+/-i*area) and N-2 ones, with spherical excess |area| mod 2pi.
%   Other polygons have no generic exact reference asserted here.
%   'operator' applies the chain of rank-two rotations in O(Nm).
%   Sphere parallel transport reference [96]:
%   https://doi.org/10.3390/sym13112092
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,H=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'finite'});
if nargin<2||isempty(points),points=eye(n,3);end
if ~isnumeric(points)||~isreal(points)||ndims(points)>2|| ...
        size(points,1)~=n||size(points,2)<3|| ...
        any(~isfinite(points(:)))|| ...
        any(abs(sqrt(sum(points.^2,1))-1)>5e-12)
    error('geodesic_holonomy:Points','POINTS must contain >=3 real unit N-vectors.');
end
m=size(points,2);
for j=1:m
    next=mod(j,m)+1;
    if 1+points(:,j)'*points(:,next)<1e-12
        error('geodesic_holonomy:Antipodal','Adjacent vertices must not be antipodal.');
    end
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('geodesic_holonomy:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    H=unitary_spectral_reference('geodesic_holonomy',n);
    H.reference='Spherical parallel transport; doi:10.3390/sym13112092';
    H.notes='Ambient product fixes the base point; restricted map is tangent holonomy.';
    if m==3
        P=points(:,1:3);
        gram=P'*P;
        volume=sqrt(max(0,det(gram)));
        denom=1+gram(1,2)+gram(2,3)+gram(3,1);
        a=2*atan2(volume,denom);
        z=[exp(1i*a);exp(-1i*a);ones(n-2,1)];
        H.finite.kind='exact';H.finite.eigenvalues=@()z;
        H.finite.moments=@(k)(n-2+2*cos(k*a))/n;
        H.notes='Exact spherical-triangle area law; orientation does not affect trace spectrum.';
    end
elseif any(strcmp(output,{'operator','implicit'}))
    H=@(X,flag)um_holonomy_apply(X,flag,points);
elseif any(strcmp(output,{'matrix','explicit'}))
    H=eye(n);
    for j=1:m
        q=points(:,mod(j,m)+1);
        p=points(:,j);
        A=q*p'-p*q';
        R=eye(n)+A+(A*A)/(1+p'*q);
        H=R*H;
    end
else
    error('geodesic_holonomy:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
