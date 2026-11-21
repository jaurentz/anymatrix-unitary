function G=um_spherical_polygon_geometry(vertices)
%UM_SPHERICAL_POLYGON_GEOMETRY Validate and triangulate geodesic polygon.
%   Vertices are M-by-3 UNIT Cartesian points on S^2 in consecutive
%   boundary order. Edges are SHORTER great-circle arcs.
%   The polygon must be simple and contained in a single well-conditioned
%   open hemisphere (all vertices dot the centroid direction > .05).
%   Concave polygons ARE supported. Holes, disconnected sets, winding
%   polygons, near-antipodal edges and vertex degeneracies are NOT.
%
%   Gnomonic projection to a tangent plane maps geodesic arcs to
%   straight line segments. An ear-clipping triangulation in that
%   plane therefore corresponds exactly to a triangulation into
%   great-circle spherical triangles. Input clockwise and anticlockwise
%   orders give the SAME region (clockwise is reversed).
%   Triangle areas are computed separately by the spherical-excess
%   atan2 formula; this is INDEPENDENT of numerical cubature.
%   [173] https://doi.org/10.1016/j.amc.2025.129335
if ~(isnumeric(vertices)&&isreal(vertices)&&ismatrix(vertices) ...
        &&size(vertices,2)==3&&size(vertices,1)>=3 ...
        &&size(vertices,1)<=24&&all(isfinite(vertices(:))))
    error('um_spherical_polygon_geometry:Vertices', ...
        'Vertices must be finite real M-by-3, 3<=M<=24.');
end
V=double(vertices);
m=size(V,1);
nr=sqrt(sum(V.^2,2));
if max(abs(nr-1))>1e-10
    error('um_spherical_polygon_geometry:Sphere', ...
        'Vertices must have unit Euclidean norm.');
end
V=V./nr;
for j=1:m
    k=mod(j,m)+1;
    if norm(V(j,:)-V(k,:))<1e-8
        error('um_spherical_polygon_geometry:Repeated', ...
            'Consecutive repeated or near-repeated vertices.');
    end
end
c=sum(V,1)';
if norm(c)<1e-8
    error('um_spherical_polygon_geometry:Hemisphere', ...
        'Centroid direction is numerically undefined.');
end
c=c/norm(c);
dots=V*c;
if min(dots)<=.05
    error('um_spherical_polygon_geometry:Hemisphere', ...
        'All vertices must lie in a well-conditioned open hemisphere.');
end
[~,axis]=min(abs(c));
reference=zeros(3,1);reference(axis)=1;
e1=reference-c*(reference'*c);e1=e1/norm(e1);
e2=cross(c,e1);
P=[(V*e1)./dots,(V*e2)./dots];
scale=max(1,max(abs(P(:))));
tol=1e-12*scale^2;
twiceArea=sum(P(:,1).*P([2:end 1],2) ...
             -P([2:end 1],1).*P(:,2));
if abs(twiceArea)<=tol
    error('um_spherical_polygon_geometry:Area', ...
        'Degenerate polygon with zero planar oriented area.');
end
if twiceArea<0
    V=flipud(V);P=flipud(P);
end
% Explicitly forbid self-intersection (unlike a winding-number mask).
for i=1:m
    j=mod(i,m)+1;
    for k=i+1:m
        l=mod(k,m)+1;
        if k==i||k==j||l==i,continue;end
        if intersection(P(i,:),P(j,:),P(k,:),P(l,:),tol)
            error('um_spherical_polygon_geometry:Intersection', ...
                'Polygon edges cross or touch away from adjacent vertices.');
        end
    end
end
% Nonadjacent repeated points are also invalid.
for i=1:m
    for j=i+1:m
        if norm(P(i,:)-P(j,:))<1e-9*scale
            error('um_spherical_polygon_geometry:Repeated', ...
                'Repeated vertices must not occur.');
        end
    end
end
% Reject immediately adjacent collinear triples: they create an
% ambiguous ear and are better represented without the middle vertex.
for i=1:m
    a=P(mod(i-2,m)+1,:);b=P(i,:);d=P(mod(i,m)+1,:);
    if abs(orient(a,b,d))<tol
        error('um_spherical_polygon_geometry:Collinear', ...
            'Remove adjacent nearly collinear vertices.');
    end
end
active=1:m;
T=zeros(m-2,3);
q=0;
while numel(active)>3
    found=false;
    h=numel(active);
    for j=1:h
        a=active(mod(j-2,h)+1);
        b=active(j);
        d=active(mod(j,h)+1);
        if orient(P(a,:),P(b,:),P(d,:))<=tol,continue;end
        ear=true;
        other=active(~ismember(active,[a b d]));
        for k=other
            x=P(k,:);
            if orient(P(a,:),P(b,:),x)>=-tol ...
                &&orient(P(b,:),P(d,:),x)>=-tol ...
                &&orient(P(d,:),P(a,:),x)>=-tol
                ear=false;break;
            end
        end
        if ear
            q=q+1;T(q,:)=[a b d];
            active(j)=[];
            found=true;break;
        end
    end
    if ~found
        error('um_spherical_polygon_geometry:Triangulation', ...
            'No valid geodesic ear found (degenerate polygon).');
    end
end
q=q+1;T(q,:)=active;
if q~=m-2
    error('um_spherical_polygon_geometry:Internal', ...
        'Unexpected number of spherical triangles.');
end
areas=zeros(m-2,1);
for j=1:m-2
    a=V(T(j,1),:)';b=V(T(j,2),:)';d=V(T(j,3),:)';
    det3=dot(a,cross(b,d));
    if det3<=1e-12
        error('um_spherical_polygon_geometry:Triangulation', ...
            'A spherical triangle is degenerate or has wrong orientation.');
    end
    areas(j)=2*atan2(det3,1+dot(a,b)+dot(b,d)+dot(d,a));
end
G=struct('vertices',V,'triangles',T,'triangle_areas',areas, ...
    'area',sum(areas),'gnomonic_vertices',P, ...
    'hemisphere_center',c,'gnomonic_frame',[e1 e2], ...
    'reversed_input',twiceArea<0, ...
    'notes',['Geodesic edges, simple polygon, possibly concave, ', ...
    'inside a single well-conditioned hemisphere; no holes. ', ...
    'Area uses spherical excess, not Gaussian cubature.']);
end

function v=orient(a,b,c)
v=(b(1)-a(1))*(c(2)-a(2)) ...
 -(b(2)-a(2))*(c(1)-a(1));
end

function yes=intersection(a,b,c,d,tol)
u=orient(a,b,c);v=orient(a,b,d);
w=orient(c,d,a);z=orient(c,d,b);
yes=(u>tol&&v<-tol||u<-tol&&v>tol) ...
    &&(w>tol&&z<-tol||w<-tol&&z>tol);
if yes,return;end
if abs(u)<=tol&&onsegment(a,b,c,tol),yes=true;return;end
if abs(v)<=tol&&onsegment(a,b,d,tol),yes=true;return;end
if abs(w)<=tol&&onsegment(c,d,a,tol),yes=true;return;end
if abs(z)<=tol&&onsegment(c,d,b,tol),yes=true;return;end
end

function yes=onsegment(a,b,p,tol)
yes=all(p>=min(a,b)-tol)&&all(p<=max(a,b)+tol);
end
