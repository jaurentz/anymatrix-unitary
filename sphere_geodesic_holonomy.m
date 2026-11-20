function [U,properties]=sphere_geodesic_holonomy(n,vertices,output)
%SPHERE_GEODESIC_HOLONOMY SO(N) parallel transport around a closed polygon.
%   U=SPHERE_GEODESIC_HOLONOMY(N,VERTICES) composes Levi-Civita
%   shortest-geodesic transports along a closed polygon on S^(N-1).
%   VERTICES is real N-by-M with M>=3 unit-vector columns, visited in
%   order v1->v2->...->vM->v1. The default is [e1 e2 e3].
%   Each arc has a unique shortest geodesic: adjacent pairs must not
%   be antipodal or numerically near antipodal. The returned SO(N)
%   ambient extension H=R(vM,v1)*...*R(v1,v2) fixes v1; restricted
%   to T_v1 S^(N-1), H is a genuine endomorphism (holonomy).
%
%   For a TRIANGLE, its vertices lie on a great 2-sphere. The
%   spherical excess has magnitude
%     Omega=2*atan2(sqrt(det(V'*V)),1+v1'*v2+v2'*v3+v3'*v1).
%   The nontrivial eigenvalues of H are exp(+/-1i*Omega);
%   the other N-2 eigenvalues equal 1. The normalized finite moments
%   are (N-2+2*cos(k*Omega))/N. In N=3 the oriented Omega has
%   numerator det([v1 v2 v3]); the multiset of eigenvalues is
%   orientation-independent. Longer loops on S^(N-1) can have
%   genuinely non-Abelian SO(N-1) holonomy; no generic phase law
%   is asserted for M>3.
%   OUTPUT='operator' composes rank-two arcs with O(M*N) matvecs.
%
%   P.-A. Absil, R. Mahony, R. Sepulchre (2008):
%   https://doi.org/10.1515/9781400830244
%   T. Needham, Visual Differential Geometry and Forms (2021):
%   https://doi.org/10.1515/9780691219899
properties={'square','unitary','normal','orthogonal','real', ...
            'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'finite'});
if nargin<2||isempty(vertices)
    vertices=zeros(n,3);
    vertices(1,1)=1;vertices(2,2)=1;vertices(3,3)=1;
end
if ~isnumeric(vertices)||~isreal(vertices)||ndims(vertices)~=2 ...
        ||size(vertices,1)~=n||size(vertices,2)<3 ...
        ||any(~isfinite(vertices(:)))
    error('sphere_geodesic_holonomy:Vertices', ...
        'VERTICES must be finite real N-by-M with M>=3.');
end
norms=sqrt(sum(vertices.^2,1));
if any(abs(norms-1)>1e-10)
    error('sphere_geodesic_holonomy:Unit', ...
        'Each VERTICES column must have unit Euclidean norm.');
end
vertices=bsxfun(@rdivide,vertices,norms);
m=size(vertices,2);
next=vertices(:,[2:m 1]);
c=sum(vertices.*next,1);
c=max(-1,min(1,c));
if any(1+c<=1e-12)
    error('sphere_geodesic_holonomy:Antipodal', ...
        'Each consecutive pair must define a well-conditioned shortest geodesic.');
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('sphere_geodesic_holonomy:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('sphere_geodesic_holonomy',n);
    U.reference='Absil, Mahony and Sepulchre (2008); Needham (2021)';
    if m==3
        V=vertices(:,1:3);
        G=V'*V;
        denominator=1+G(1,2)+G(2,3)+G(3,1);
        numerator=sqrt(max(0,det(G)));
        omega=2*atan2(numerator,denominator);
        U.finite.kind='exact';
        U.finite.eigenvalues=@()[exp(1i*omega);exp(-1i*omega);ones(n-2,1)];
        U.finite.moments=@(k)(n-2+2*cos(k*omega))/n;
        if abs(omega)<5e-14
            U.finite.atoms=struct('angles',0,'weights',1);
        else
            U.finite.atoms=struct('angles',[omega;-omega;0], ...
                                  'weights',[1;1;n-2]/n);
        end
        U.notes=['Triangle holonomy: exact spherical-excess conjugate ' ...
            'pair and N-2 unit eigenvalues. In higher dimensions ' ...
            'the sign of the excess depends on orientation choice, ' ...
            'but its eigenvalue multiset does not.'];
    else
        U.notes=['Closed polygon is an SO(N) ambient lift fixing its ' ...
            'basepoint, hence a genuine tangent-space holonomy. ' ...
            'Non-Abelian higher-dimensional polygon spectra may vary; ' ...
            'no general finite phase law is asserted for M>3.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('sphere_geodesic_holonomy:Mode', ...
          'Supported modes: matrix, operator, spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_spherical_holonomy_apply(X,flag,vertices,c);
else
    U=eye(n);
    for j=1:m
        x=vertices(:,j);y=next(:,j);
        K=y*x'-x*y';
        R=eye(n)+K+(K*K)/(1+c(j));
        U=R*U;
    end
end
end
