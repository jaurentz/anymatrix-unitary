function [U,properties]=mesh_spin_dirac(vertices,faces,t,mode,gauges,boundary,steps)
%MESH_SPIN_DIRAC Face-quaternionic mesh Dirac unitary propagator.
%   U=MESH_SPIN_DIRAC(V,F,T,MODE,GAUGES,BOUNDARY,STEPS) builds the
%   face-based INTRINSIC spin-Dirac quaternionic matrix of
%   Ye, Diamanti, Tang, Guibas and Hoffmann, Computer Graphics Forum 37
%   (2018), Eq. (19), DOI 10.1111/cgf.13494:
%
%      (D0 phi)_i = (1/2)*sum_(j adjacent to i) E_ij*cos(theta_ij/2)*phi_j,
%      E_ij = 2 H_ij + e_ij, H_ij=|e_ij|*tan(theta_ij/2)/2.
%
%   The per-face spinor is a quaternion q=[q0 q1 q2 q3], and
%   left multiplication by E_ij is represented by the 4-by-4
%   REAL matrix L(E_ij). On an oriented triangle mesh, E_ji=conj(E_ij)
%   and L(conj(q))=L(q)', hence D0'=D0 exactly in mathematical
%   arithmetic. This is an ACTUAL face-spinor construction, NOT
%   SIMPLICIAL_HODGE_DIRAC or the naive flat-torus two-spinor lattice.
%
%   Input V is an nV-by-3 real vertex array, F an nF-by-3
%   consistently oriented triangle array. All vertices must be used.
%   For default V,F, use a planar 2-triangle unit square:
%     V=[(0,0),(1,0),(1,1),(0,1)] and F=[1 2 3;1 3 4].
%   Each shared edge has incident faces i,j with opposite induced
%   oriented edge. Take e_ij along FACE i's oriented boundary and
%     theta_ij=atan2(dot(e_ij/|e|,cross(n_i,n_j)),dot(n_i,n_j)).
%   Theta is the SIGNED normal bending angle in (-pi,pi). With this
%   convention E_ij*cos(theta/2)
%     = |e|*sin(theta/2)+e_ij*cos(theta/2).
%   This is a fixed signed angle convention for a planar classical
%   face edge-constrained net. It differs from simply putting a
%   nonnegative absolute dihedral angle into an unoriented edge;
%   the orientation is specified so the operator is reproducible.
%
%   AREA METRIC (explicitly DERIVED, NOT Ye et al. Eq. (19)):
%   Use positive lumped face areas M=diag(repelem(area_i,4)).
%   Interpret D_phys=M^-1*D0 as a face-area generalized Dirac
%   D0*phi=lambda*M*phi. Then D_phys'*M=M*D_phys and the
%   Euclidean-normalized Hermitian H=M^-1/2*D0*M^-1/2
%   gives U(T)=expm(-1i*T*H). The physical propagation is
%   Vphys=M^-1/2*U*M^1/2 and is M-unitary, not generally Euclidean
%   unitary. Ye et al. also use face areas in a DIFFERENT squared/
%   averaged generalized eigenproblem, Eq. (27): do not identify
%   the two discretizations. This positive mass is a benchmark choice.
%
%   GAUGES optional nF-by-4 rows of UNIT quaternions g_i, in
%   [scalar i j k] order. G=diag(L(g_i)); D0_g=G'*D0*G.
%   This implements exactly the local quaternionic change of frame.
%   It is a LOCAL GAUGE of the chosen spin structure, not a method
%   for constructing the globally distinct spin structures of a torus.
%   The face metric commutes with G, so Hg=G'*H*G and Ug=G'*U*G.
%
%   BOUNDARY='free' (default): boundary edges have no adjacent face
%   contribution (natural omitted-neighbor boundary); 'closed'
%   rejects any boundary edges. The paper does not prescribe a
%   self-adjoint continuum boundary-value problem here: the free
%   convention is the explicitly stated finite discrete one.
%
%   MODE='matrix' (default): dense true exponential, 4*nF<=256.
%        'model' : dense U, sparse H, D0, D_phys, M, physical
%                  M-unitary Vphys, geometry and gauge, 4*nF<=256.
%        'generator' : sparse matrices H,D0,D_phys,M and geometry,
%                  without forming expm, 4*nF<=12000.
%        'cayley_operator' : matrix-free in the sense of NO dense U,
%                  sparse rational Cayley approximation to exp(-itH).
%                  Cached sparse LU factors; not the exact exponential.
%                  Steps default 16. Forward, adjoint, ordinary
%                  transpose and multiple RHS with size flag.
%        'spectral' : exact finite phase multiset for the default
%                  TWO adjacent triangles (4 copies of each +/-λ),
%                  or for T=0. Generic closed mesh phase law UNKNOWN.
%                  Real Dirac nodes ±λ are NOT unitary eigenphases.
%
%   Mathematically U is unitary. Numerical expm and rational solves
%   have floating point errors. A Cayley step is itself unitary in
%   exact arithmetic but a finite number of steps is only an
%   O(T^3/STEPS^2) approximation to the true exponential at fixed H.
%   Matrix-free approximation is not advertised as true expm.
%
%   [170] Ye et al. (2018), Eq. (12), (19), (27), Theorem 4.3:
%   https://doi.org/10.1111/cgf.13494
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
defaultMesh=isempty(vertices);
if defaultMesh
    vertices=[0 0 0;1 0 0;1 1 0;0 1 0];
end
if nargin<2||isempty(faces)
    if ~defaultMesh
        error('mesh_spin_dirac:Faces', ...
            'FACES must be supplied for custom vertices.');
    end
    faces=[1 2 3;1 3 4];
end
if nargin<3||isempty(t),t=.7;end
if nargin<4||isempty(mode),mode='matrix';end
if nargin<5,gauges=[];end
if nargin<6||isempty(boundary),boundary='free';end
if nargin<7||isempty(steps),steps=16;end
if ~(isnumeric(t)&&isscalar(t)&&isreal(t)&&isfinite(t))
    error('mesh_spin_dirac:Time','T must be a finite real scalar.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('mesh_spin_dirac:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','generator', ...
        'cayley_operator','spectral'}))
    error('mesh_spin_dirac:Mode','Unsupported output MODE.');
end
if ~(ischar(boundary)||(isstring(boundary)&&isscalar(boundary)))
    error('mesh_spin_dirac:Boundary','BOUNDARY must be scalar text.');
end
boundary=lower(char(boundary));
if ~any(strcmp(boundary,{'free','closed'}))
    error('mesh_spin_dirac:Boundary','Use free or closed.');
end
if ~(isnumeric(steps)&&isscalar(steps)&&isreal(steps) ...
        &&isfinite(steps)&&steps==fix(steps)&&steps>=1&&steps<=4096)
    error('mesh_spin_dirac:Steps','STEPS must be integer 1..4096.');
end
if ~(isnumeric(vertices)&&isreal(vertices)&&ismatrix(vertices) ...
        &&size(vertices,2)==3 && size(vertices,1)>=3 ...
        &&all(isfinite(vertices(:))))
    error('mesh_spin_dirac:Vertices', ...
        'V must be a finite real nV-by-3 vertex array.');
end
V=double(vertices);nv=size(V,1);
if ~(isnumeric(faces)&&isreal(faces)&&ismatrix(faces) ...
        &&size(faces,2)==3 && size(faces,1)>=1 ...
        &&all(isfinite(faces(:)))&&all(faces(:)>=1) ...
        &&all(faces(:)<=nv)&&all(faces(:)==fix(faces(:))))
    error('mesh_spin_dirac:Faces', ...
        'F must be an nF-by-3 integer array of vertex indices.');
end
F=double(faces);nf=size(F,1);N=4*nf;
if N>12000
    error('mesh_spin_dirac:Order','Spinor dimension must be <=12000.');
end
if any(F(:,1)==F(:,2)|F(:,2)==F(:,3)|F(:,3)==F(:,1)) ...
        ||size(unique(sort(F,2),'rows'),1)~=nf ...
        ||numel(unique(F(:)))~=nv
    error('mesh_spin_dirac:Topology', ...
        'All vertices must be used and faces must be distinct triangles.');
end
p0=V(F(:,1),:);p1=V(F(:,2),:);p2=V(F(:,3),:);
normals=cross(p1-p0,p2-p0,2);
area=.5*sqrt(sum(normals.^2,2));
span=max(1,max(max(V)-min(V)));
if any(~isfinite(area)|area<=1e-13*span^2)
    error('mesh_spin_dirac:Geometry','Degenerate triangle area.');
end
unitNormals=normals./(2*area);
faceEdges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
faceIds=repmat((1:nf)',3,1);
half=[sort(faceEdges,2),faceEdges,faceIds];
half=sortrows(half,[1 2]);
D0=spalloc(N,N,32*nf);
boundaryCount=0;
edgeCount=0;
edgeI=zeros(3*nf,1);edgeJ=zeros(3*nf,1);
edgeLengths=zeros(3*nf,1);edgeAngles=zeros(3*nf,1);
k=1;
while k<=size(half,1)
    j=k;
    while j<size(half,1)&&all(half(j+1,1:2)==half(k,1:2))
        j=j+1;
    end
    if j-k>=2
        error('mesh_spin_dirac:Manifold', ...
            'An edge may belong to at most two faces.');
    end
    if j==k
        boundaryCount=boundaryCount+1;
    else
        if half(k,3)~=half(j,4) || half(k,4)~=half(j,3)
            error('mesh_spin_dirac:Orientation', ...
                'Adjacent faces must induce opposite edge directions.');
        end
        i=half(k,5);l=half(j,5);
        v=V(half(k,4),:)-V(half(k,3),:);
        lengthEdge=norm(v);
        if lengthEdge<=1e-13*span
            error('mesh_spin_dirac:Geometry','Zero length edge.');
        end
        ci=unitNormals(i,:);cj=unitNormals(l,:);
        sine=dot(v/lengthEdge,cross(ci,cj));
        cosine=max(-1,min(1,dot(ci,cj)));
        theta=atan2(sine,cosine);
        if abs(cos(theta/2))<=1e-12
            error('mesh_spin_dirac:Fold', ...
                'Opposite face normals are unsupported.');
        end
        % Source Eq. (12) & (19), with explicit signed edge convention.
        quaternion=[lengthEdge*sin(theta/2),v*cos(theta/2)];
        B=.5*quatleft(quaternion);
        I=4*(i-1)+(1:4);J=4*(l-1)+(1:4);
        D0(I,J)=D0(I,J)+B;
        D0(J,I)=D0(J,I)+B.';
        edgeCount=edgeCount+1;
        edgeI(edgeCount)=i;edgeJ(edgeCount)=l;
        edgeLengths(edgeCount)=lengthEdge;
        edgeAngles(edgeCount)=theta;
    end
    k=j+1;
end
if strcmp(boundary,'closed') && boundaryCount>0
    error('mesh_spin_dirac:Boundary', ...
        'Closed boundary convention forbids boundary edges.');
end
if isempty(gauges)
    gauges=repmat([1 0 0 0],nf,1);
end
if ~(isnumeric(gauges)&&isreal(gauges)&&ismatrix(gauges) ...
        &&isequal(size(gauges),[nf 4]) ...
        &&all(isfinite(gauges(:))))
    error('mesh_spin_dirac:Gauges', ...
        'GAUGES must be an nF-by-4 real array of unit quaternions.');
end
if max(abs(sum(gauges.^2,2)-1))>5e-12
    error('mesh_spin_dirac:Gauges', ...
        'Every gauge quaternion must be unit length.');
end
G=spalloc(N,N,16*nf);
for i=1:nf
    I=4*(i-1)+(1:4);
    G(I,I)=quatleft(gauges(i,:));
end
D0=G.'*D0*G;
mass=repelem(area,4);
M=spdiags(mass,0,N,N);
Ri=spdiags(1./sqrt(mass),0,N,N);
H=Ri*D0*Ri;
Dphys=spdiags(1./mass,0,N,N)*D0;
edges=[edgeI(1:edgeCount),edgeJ(1:edgeCount)];
common=struct('H',H,'D_integrated',D0,'D_physical',Dphys, ...
    'metric',M,'face_areas',area,'face_normals',unitNormals, ...
    'edge_face_pairs',edges, ...
    'edge_lengths',edgeLengths(1:edgeCount), ...
    'edge_signed_angles',edgeAngles(1:edgeCount), ...
    'boundary_edge_count',boundaryCount,'boundary',boundary, ...
    'gauges',gauges,'gauge_matrix',G, ...
    'faces',F,'vertices',V,'spinor_components_per_face',4, ...
    'dimension',N,'time',t, ...
    'notes',['Ye et al. intrinsic quaternionic FACE Dirac Eq (19), ', ...
      'four REAL spinor coordinates per face. The per-face AREA ', ...
      'generalized mass is an explicitly derived separate model; ', ...
      'not the author''s Eq (27) squared/averaged eigenproblem. ', ...
      'Gauge is local quaternion frame, not global spin topology. ', ...
      'Boundary missing-neighbor convention is a finite choice. ', ...
      'Generic finite Dirac eigenvalues and zero modes unknown.']);
if strcmp(mode,'generator'),U=common;return;end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('mesh_spin_dirac',N);
    U.reference='Ye et al., Comput. Graph. Forum 37 (2018), Eq. (12),(19)';
    U.boundary=boundary;U.time=t;
    if t==0
        z=ones(N,1);
    elseif nf==2 && edgeCount==1
        % Source blocks [0,B;B',0] have singular values |e|/2,
        % each repeated four times. Metric normalization divides
        % by sqrt(area1*area2) and gauge cannot change eigenvalues.
        lambda=edgeLengths(1)/(2*sqrt(prod(area)));
        z=[exp(-1i*t*lambda)*ones(4,1); ...
           exp( 1i*t*lambda)*ones(4,1)];
    else
        z=[];
    end
    if ~isempty(z)
        U.finite.kind='exact';
        U.finite.eigenvalues=@()z;
        U.finite.moments=@(k)moments(z,k);
        U.finite.atoms=struct('angles',angle(z), ...
            'weights',ones(N,1)/N);
        U.notes=['Exact two-face bipartite quaternion edge Dirac ', ...
            'phase multiplicity four for +/-lambda, or time-zero. ', ...
            'This is the AREA-normalized derived metric model.'];
    else
        U.notes=['Generic immersed triangle-mesh spin-Dirac ', ...
            'finite eigenphases are UNKNOWN. Quadrature geometric ', ...
            'data are not unitary eigenphases.'];
    end
    return
end
if strcmp(mode,'cayley_operator')
    h=t/steps;
    plusMatrix=speye(N)+.5i*h*H;
    minusMatrix=speye(N)-.5i*h*H;
    plus=decomposition(plusMatrix,'lu');
    minus=decomposition(minusMatrix,'lu');
    U=@(X,flag)um_mesh_spin_dirac_cayley(X,flag, ...
        N,steps,plus,minus,plusMatrix,minusMatrix);
    return
end
if N>256
    error('mesh_spin_dirac:ExplicitOrder', ...
        'Dense matrix/model expm requires 4*nF<=256.');
end
Q=expm(-1i*t*full(H));
if strcmp(mode,'matrix'),U=Q;return;end
R=spdiags(sqrt(mass),0,N,N);
common.U=Q;
common.V_physical=Ri*Q*R;
U=common;
end

function L=quatleft(q)
% Left quaternion multiplication for [1,i,j,k] column coordinates.
a=q(1);b=q(2);c=q(3);d=q(4);
L=[a -b -c -d;b a -d c;c d a -b;d -c b a];
end

function m=moments(z,k)
if ~(isnumeric(k)&&isreal(k)&&all(isfinite(k(:))) ...
        &&all(k(:)==fix(k(:))))
    error('mesh_spin_dirac:Moments','Moment orders must be integers.');
end
m=reshape(arrayfun(@(r)mean(z.^r),k(:)),size(k));
end
