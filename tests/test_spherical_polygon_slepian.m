function tests=test_spherical_polygon_slepian
%TEST_SPHERICAL_POLYGON_SLEPIAN Geodesic polygon positive cubature.
tests=functiontests(localfunctions);
end

function test_triangle_square_concave_reference_and_spectrum(tc)
shapes={[-.45 -.32; .64 -.25; .37 .53], ...
        [-.53 -.43;.52 -.43;.52 .45;-.53 .45], ...
        [-.75 -.48;.72 -.46;.09 -.11;.72 .52;-.72 .49]};
for h=1:numel(shapes)
    V=unit_north(shapes{h});
    for L=[1 2 4 6]
        n=L^2;
        P=spherical_polygon_slepian(L,V,1e-9,'parameters');
        K=P.concentration;Q=P.basis;
        lam=P.concentration_eigenvalues;
        verifyEqual(tc,size(Q),[n n]);
        verifyTrue(tc,isreal(Q));
        verifyTrue(tc,P.quadrature_converged);
        verifyFalse(tc,P.quadrature_certified);
        verifyGreaterThan(tc,P.quadrature_levels,1);
        verifyLessThan(tc,P.quadrature_difference_fro, ...
            1e-8*max(1,P.shannon_number));
        verifyLessThan(tc,P.independent_area_difference, ...
            1e-8*max(1,P.geometric_area));
        verifyLessThan(tc,norm(K-K','fro'),3e-12*n);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),3e-10*n);
        verifyLessThan(tc,norm(Q'*K*Q-diag(lam),'fro'), ...
            3e-9*n);
        verifyLessThan(tc,abs(trace(K)-P.shannon_number), ...
            5e-9*n);
        verifyGreaterThanOrEqual(tc,min(lam),-3e-10);
        verifyLessThanOrEqual(tc,max(lam),1+3e-10);
        verifyLessThan(tc,norm(spherical_polygon_slepian(L,V) ...
            -Q','fro'),5e-11*n);
        if L<=4
            % Independent gnomonic planar triangle integral of
            % physically evaluated Condon--Shortley harmonics.
            R=independent_gnomonic_quadrature(L,V,36);
            verifyLessThan(tc,norm(K-R,'fro'), ...
                7e-10*n);
        end
        if L==1
            verifyLessThan(tc,abs(K-P.geometric_area/(4*pi)), ...
                2e-11);
        end
    end
end
end

function test_geodesic_area_from_independent_spherical_excess(tc)
V=unit_north([-.6 -.4;.58 -.4;.13 -.12;.65 .57;-.6 .5]);
G=um_spherical_polygon_geometry(V);
verifyEqual(tc,size(G.triangles,1),size(V,1)-2);
sumArea=0;
for k=1:size(G.triangles,1)
    Vtri=G.vertices(G.triangles(k,:),:);
    side=[acos(dot(Vtri(2,:),Vtri(3,:))), ...
          acos(dot(Vtri(1,:),Vtri(3,:))), ...
          acos(dot(Vtri(1,:),Vtri(2,:)))];
    s=sum(side)/2;
    val=prod(tan([s,s-side]/2));
    area=4*atan(sqrt(max(0,val))); % independent L'Huilier formula
    sumArea=sumArea+area;
end
verifyLessThan(tc,abs(sumArea-G.area),2e-11);
verifyGreaterThan(tc,G.area,.1);
verifyLessThan(tc,G.area,2*pi);
end

function test_orientation_rotation_and_geodesic_covariance(tc)
V=unit_north([-.48 -.34;.51 -.27;.65 .19;.2 .58;-.62 .43]);
angles=[.31 .7 -.22];
R=physical_rotation(angles);
Vr=(R*V.').';
for L=[2 4 6]
    A=spherical_polygon_slepian(L,V,2e-9,'parameters');
    B=spherical_polygon_slepian(L,flipud(V),2e-9,'parameters');
    C=spherical_polygon_slepian(L,Vr,2e-9,'parameters');
    verifyLessThan(tc,norm(A.concentration-B.concentration,'fro'), ...
        3e-10*L^2);
    verifyLessThan(tc,abs(A.geometric_area-B.geometric_area), ...
        2e-13);
    verifyLessThan(tc,abs(A.geometric_area-C.geometric_area), ...
        5e-13);
    verifyLessThan(tc,norm(sort(A.concentration_eigenvalues) ...
        -sort(C.concentration_eigenvalues)), ...
        5e-9*L^2);
end
end

function test_geometry_and_unknown_eigenphases(tc)
V=unit_north([-.45 -.32;.64 -.25;.37 .53]);
G=spherical_polygon_slepian(3,V,[], 'geometry');
S=spherical_polygon_slepian(3,V,[], 'spectral');
verifyEqual(tc,G.triangles,[1 2 3]);
verifyGreaterThan(tc,G.area,0);
verifyEqual(tc,S.dimension,9);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
end

function test_validation_reject_ambiguous_domains_and_uncertainty(tc)
[U,props]=spherical_polygon_slepian();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(props,'orthogonal')));
V=unit_north([-.5 -.35;.6 -.3;.5 .45;-.5 .4]);
verifyError(tc,@()spherical_polygon_slepian(11,V), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()spherical_polygon_slepian(3,V,-1), ...
    'spherical_polygon_slepian:Tolerance');
verifyError(tc,@()spherical_polygon_slepian(3,V,1e-9,'operator'), ...
    'spherical_polygon_slepian:Mode');
verifyError(tc,@()spherical_polygon_slepian(3,V,1e-25), ...
    'spherical_polygon_slepian:Quadrature');
bad=V;bad(2,:)=bad(1,:);
verifyError(tc,@()spherical_polygon_slepian(3,bad), ...
    'um_spherical_polygon_geometry:Repeated');
bad=V;bad(2,:)=2*bad(2,:);
verifyError(tc,@()spherical_polygon_slepian(3,bad), ...
    'um_spherical_polygon_geometry:Sphere');
% Self-intersection in projected XY plane, with nonzero signed area.
bow=unit_north([-.71 -.62;.64 .67;-.53 .58;.67 -.46]);
verifyError(tc,@()spherical_polygon_slepian(3,bow), ...
    'um_spherical_polygon_geometry:Intersection');
wide=[0 0 1;1 0 0;0 0 -1;-1 0 0];
verifyError(tc,@()spherical_polygon_slepian(3,wide), ...
    'um_spherical_polygon_geometry:Hemisphere');
% Adjacent geodesic edge collinearity is intentionally unsupported.
flat=unit_north([-.5 -.4;0 -.4;.5 -.4;.5 .5;-.5 .5]);
verifyError(tc,@()spherical_polygon_slepian(3,flat), ...
    'um_spherical_polygon_geometry:Collinear');
end

function V=unit_north(P)
V=[P,ones(size(P,1),1)];
V=V./sqrt(sum(V.^2,2));
end

function R=physical_rotation(a)
x=a(1);y=a(2);z=a(3);
C=[cos(x),-sin(x),0;sin(x),cos(x),0;0,0,1];
Y=[cos(y),0,sin(y);0,1,0;-sin(y),0,cos(y)];
D=[cos(z),-sin(z),0;sin(z),cos(z),0;0,0,1];
R=C*Y*D;
end

function K=independent_gnomonic_quadrature(L,V,q)
% This alternative calculation uses only the gnomonic x/z,y/z plane
% Jacobian (1+x²+y²)^(-3/2), not the spherical chord
% det[a,b,c]/||x||³ Jacobian or the production Gaussian rule.
G=um_spherical_polygon_geometry(V);
j=(1:q-1)';
off=j./sqrt(4*j.^2-1);
J=diag(off,1)+diag(off,-1);
[T,D]=eig(J);
[x,ind]=sort(real(diag(D)));
u=(x+1)/2;
w=real(T(1,ind).').^2;
S=repmat(u,1,q);H=repmat(u.',q,1);
U=S(:);H=H(:);
wp=repmat(w,1,q).*repmat(w.',q,1).*S;
n=L^2;K=zeros(n);
for tri=1:size(G.triangles,1)
    a=G.vertices(G.triangles(tri,1),:);
    b=G.vertices(G.triangles(tri,2),:);
    c=G.vertices(G.triangles(tri,3),:);
    P=[a(1:2)/a(3);b(1:2)/b(3);c(1:2)/c(3)];
    xyr=P(1,:).*(1-U)+P(2,:).*(U.*(1-H)) ...
        +P(3,:).*(U.*H);
    xy=[xyr(:,1),xyr(:,2)];
    det2=abs(det([P(2,:)-P(1,:);P(3,:)-P(1,:)]));
    norm3=sqrt(1+sum(xy.^2,2));
    coords=[xy,ones(size(xy,1),1)]./norm3;
    weights=wp(:)*det2./norm3.^3;
    B=zeros(q^2,n);
    for node=1:q^2
        v=coords(node,:);
        phi=atan2(v(2),v(1));
        idx=0;
        for m=0:L-1
            if m==0,kinds=0;else,kinds=[1 2];end
            for kind=kinds
                for ell=m:L-1
                    Pm=legendre(ell,v(3));
                    f=sqrt((2*ell+1)/2*exp( ...
                        gammaln(ell-m+1)-gammaln(ell+m+1)));
                    if kind==0,angular=1/sqrt(2*pi);
                    elseif kind==1,angular=cos(m*phi)/sqrt(pi);
                    else,angular=sin(m*phi)/sqrt(pi);end
                    idx=idx+1;B(node,idx)=f*Pm(m+1)*angular;
                end
            end
        end
    end
    K=K+B'*(weights.*B);
end
end
