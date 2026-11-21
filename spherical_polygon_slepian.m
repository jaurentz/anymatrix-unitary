function [U,properties]=spherical_polygon_slepian(L,vertices,tol,mode)
%SPHERICAL_POLYGON_SLEPIAN Real orthogonal Slepian basis on a geodesic polygon.
%   U=SPHERICAL_POLYGON_SLEPIAN(L,VERTICES,TOL,MODE) constructs the
%   N-by-N ANALYSIS matrix Q', N=L^2, of the bandlimited spherical
%   concentration K_ab=int_polygon Y_a(Omega)Y_b(Omega) dOmega.
%   K is REAL SYMMETRIC and a positive contraction in exact arithmetic,
%   but it is NOT orthogonal/unitary; its concentration eigenvalues
%   are NOT the eigenphases of the orthogonal matrix Q'.
%
%   VERTICES is an M-by-3 array of unit Cartesian vectors defining a
%   simple geodesic polygon with SHORTER great-circle edges.
%   Clockwise/counterclockwise and convex/concave polygons supported
%   (but all vertices must lie in the same well-conditioned open
%   hemisphere). Holes, disconnected regions, winding/self-crossing
%   edges, antipodal edges and patches wider than a hemisphere are
%   rejected rather than silently assigned ambiguous domains.
%   At most 24 vertices and L<=10 for this dense reference.
%
%   Gnomonic projection sends geodesic arcs to planar line segments.
%   Ear clipping triangulates even CONCAVE simple polygons into
%   spherical geodesic triangles. On each triangle, radial projection
%   of a planar affine barycentric patch has exact surface Jacobian:
%     x(s,t)=(1-s)a+s(1-t)b+st c,
%     dOmega=s*det([a,b,c])/norm(x)^3 ds dt.
%   Gaussian tensor quadrature computes K. No finite q gives
%   an EXACT integral in general: K approximates the mathematical
%   concentration, and Q' approximates its eigenbasis.
%
%   TOL (default 1e-9): requested absolute FROBENIUS/NORMALIZED
%   convergence proxy on K, obtained from successive Gaussian orders,
%   and separate area/spherical-addition-trace checks. Increase
%   quadrature order until these diagnostics are <=TOL scale or
%   throw spherical_polygon_slepian:Quadrature. The difference of
%   successive rules is an a posteriori CONVERGENCE INDICATOR,
%   NOT a rigorous upper error bound on the true integral. No
%   certified absolute spectral bound or ARPIST verification claimed.
%
%   MODE='matrix' or 'analysis': Q' (default, real orthogonal)
%        'basis': Q (orthogonal)
%        'concentration': K (NONUNITARY, approximate integral)
%        'parameters': K,Q,lambda,exact geometric area, integration
%                      diagnostics, basis labels and tie tolerance
%        'geometry': triangulation, exact spherical-excess area,
%                    gnomonic chart diagnostics (no quadrature)
%        'spectral': generic finite analysis eigenphases UNKNOWN.
%
%   Harmonics are real Condon--Shortley with the SAME m=0, then
%   cos(m*phi)/sin(m*phi), ell=m,...,L-1 basis as spherical_slepian.
%   Trace(K) should converge to N*area/(4*pi) (addition theorem).
%   Positive numerical cubature yields PSD K up to floating roundoff.
%   The returned basis is numerically orthogonal after eig(K);
%   it is NOT claimed to be an exact basis of the mathematical K
%   before quadrature convergence. Nearly repeated eigenspaces use
%   deterministic pivot projection with explicit tie tolerance.
%
%   Dense O(N^3+M*q^2*N^2) setup and O(N^2) storage/action;
%   q<=64, M<=24. No matrix-free accelerated polygon transform.
%
%   [156] https://doi.org/10.1137/S0036144504445765
%   [172] https://doi.org/10.1016/j.cam.2022.114822 (radial projection)
%   [173] https://doi.org/10.1016/j.amc.2025.129335 (spherical polygons)
properties={'square','unitary','orthogonal','real','normal', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','real','positive','finite','<=',10});
if nargin<2||isempty(vertices)
    vertices=[0 0 1;.6 0 .8;.3 .6 .8];
    vertices=vertices./sqrt(sum(vertices.^2,2));
end
if nargin<3||isempty(tol),tol=1e-9;end
if ~(isnumeric(tol)&&isscalar(tol)&&isreal(tol)&&isfinite(tol) ...
        &&tol>0&&tol<=.1)
    error('spherical_polygon_slepian:Tolerance', ...
        'TOL must be finite real, 0<TOL<=0.1.');
end
if nargin<4||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('spherical_polygon_slepian:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','analysis','basis', ...
        'concentration','parameters','spectral','geometry'}))
    error('spherical_polygon_slepian:Mode','Unsupported MODE.');
end
G=um_spherical_polygon_geometry(vertices);
n=L^2;
if strcmp(mode,'geometry'),U=G;return;end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('spherical_polygon_slepian',n);
    U.geometric_area=G.area;
    U.shannon_number=n*G.area/(4*pi);
    U.polygon_vertices=G.vertices;
    U.reference='Li--Jiao 2023; Sommariva--Vianello 2025';
    U.notes=['Finite real orthogonal eigenbasis of approximated ', ...
        'spherical polygon concentration. Generic eigenphases ', ...
        'UNKNOWN; concentration factors lambda are NOT phases ', ...
        'of the orthogonal Q transpose.'];
    return
end
targetTrace=n*G.area/(4*pi);
q=max(8,L+3);
[Kold,areaOld]=um_spherical_polygon_quadrature(L,G,q);
converged=false;
indicator=Inf;traceError=Inf;areaError=Inf;
levels=1;
while 2*q<=64
    q=2*q;
    [K,areaNumerical]=um_spherical_polygon_quadrature(L,G,q);
    levels=levels+1;
    indicator=norm(K-Kold,'fro');
    traceError=abs(trace(K)-targetTrace);
    areaError=abs(areaNumerical-G.area);
    scale=max(1,targetTrace);
    if max([indicator,traceError,areaError])<=tol*scale
        converged=true;break;
    end
    Kold=K;areaOld=areaNumerical;
end
if ~converged
    error('spherical_polygon_slepian:Quadrature', ...
        ['Gaussian polygon cubature did not satisfy TOL=%g ', ...
         'by q=%d; differences/trace/area=%g/%g/%g. ', ...
         'These are empirical indicators, not certified bounds.'], ...
        tol,q,indicator,traceError,areaError);
end
symDefect=norm(K-K','fro');
if symDefect>tol*max(1,norm(K,'fro'))
    error('spherical_polygon_slepian:Symmetry', ...
        'Quadrature concentration symmetry defect exceeds tolerance.');
end
K=(K+K')/2; % only after symmetry defect was measured and checked
[V,D]=eig(K);
[lambda,perm]=sort(real(diag(D)),'descend');
if min(lambda)<-max(1e-10,20*tol) ...
        ||max(lambda)>1+max(1e-10,20*tol)
    error('spherical_polygon_slepian:Contraction', ...
        'Quadrature K falls outside positive-contraction tolerance.');
end
Q=real(V(:,perm));
tieTolerance=2e-12*max(1,norm(K,2));
i=1;
while i<=n
    j=i;
    while j<n&&abs(lambda(j)-lambda(j+1))<=tieTolerance
        j=j+1;
    end
    if j>i
        % Complete tie subspace is already orthonormal from eig(K).
        % Pivoted QR of its coordinate restriction selects well-
        % conditioned canonical axes, unlike sequential Gram-Schmidt
        % of full projector columns (which can amplify roundoff by
        % 1e-6 for tiny clustered Slepian eigenvalues).
        Z=Q(:,i:j);
        [~,~,pivots]=qr(Z','vector');
        rankCount=j-i+1;
        W=Z*(Z(pivots(1:rankCount),:))';
        [W,~]=qr(W,0);
        Q(:,i:j)=W;
    end
    i=j+1;
end
for k=1:n
    [~,pivot]=max(abs(Q(:,k)));
    if Q(pivot,k)<0,Q(:,k)=-Q(:,k);end
end
degree=zeros(n,1);order=zeros(n,1);kind=zeros(n,1);h=0;
for m=0:L-1
    if m==0,types=0;else,types=[1 2];end
    for k=types
        for ell=m:L-1
            h=h+1;degree(h)=ell;order(h)=m;kind(h)=k;
        end
    end
end
switch mode
    case {'matrix','analysis'}
        U=Q.';
    case 'basis'
        U=Q;
    case 'concentration'
        U=K;
    case 'parameters'
        U=struct('analysis',Q.','basis',Q,'concentration',K, ...
            'concentration_eigenvalues',lambda, ...
            'bandlimit',L,'vertices',G.vertices, ...
            'triangles',G.triangles,'geometric_area',G.area, ...
            'numerical_area',areaNumerical, ...
            'shannon_number',targetTrace, ...
            'quadrature_order',q,'quadrature_levels',levels, ...
            'quadrature_difference_fro',indicator, ...
            'independent_area_difference',areaError, ...
            'addition_trace_difference',traceError, ...
            'quadrature_requested_tolerance',tol, ...
            'quadrature_converged',converged, ...
            'quadrature_certified',false, ...
            'symmetry_defect_before_averaging',symDefect, ...
            'degree',degree,'azimuthal_order',order, ...
            'harmonic_kind',kind,'eigenvalue_original_index',perm(:), ...
            'eigenvalue_tie_tolerance',tieTolerance, ...
            'notes',['Positive adaptive-in-order Gaussian/Duffy ', ...
            'cubature on gnomonic-ear spherical polygon triangles; ', ...
            'successive-rule differences are NOT rigorous error ', ...
            'bounds. Mathematical K is a concentration contraction ', ...
            'and NOT unitary. Q transpose is numerically orthogonal.']);
end
end
