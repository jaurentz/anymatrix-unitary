function [U,properties]=simplicial_hodge_dirac(vertices,faces,t,mode,steps)
%SIMPLICIAL_HODGE_DIRAC Metric-normalized discrete Hodge–Dirac propagator.
%   U=SIMPLICIAL_HODGE_DIRAC(V,F,T) returns the mathematically EXACT
%   finite Euclidean unitary exp(-1i*T*H), from oriented triangles
%   F and embedded Cartesian vertices V, via the REAL SYMMETRIC
%   Hodge–de Rham / KÄHLER–DIRAC operator on simplicial cochains.
%   IMPORTANT: this is NOT a 2-component spinor Dirac operator,
%   and NOT Ye–Diamanti–Tang–Guibas–Hoffmann's intrinsic/extrinsic
%   spin-connection mesh Dirac discretization. It is a distinct,
%   topological discrete exterior-calculus construction adjacent
%   to open spin-Dirac issue #59, which it does NOT resolve.
%
%   Cochain degrees: C^0 of V vertices, C^1 of E oriented edges
%   (always ascending vertex index), C^2 of F oriented triangles.
%   d0:C^0->C^1 and d1:C^1->C^2 have d1*d0=0 EXACTLY
%   in integer arithmetic. Let d combine these two coboundaries.
%   Use explicit POSITIVE diagonal cochain inner products
%     M0(v) = sum_{f incident on v} area(f)/3,
%     M1(e) = sum_{f incident on e} area(f)/(3*length(e)^2),
%     M2(f) = 1/area(f).
%   These are a documented simple LUMPED POSITIVE metric, not
%   the exact Whitney mass matrix or a proven convergent DEC
%   circumcentric Hodge star. On curved meshes their geometric
%   approximation quality is a separate numerical question.
%
%   Physical D=d+M^(-1)*d'*M is M-selfadjoint (D'*M=M*D).
%   Euclidean H=M^(1/2)*D*M^(-1/2)
%              =M^(1/2)*d*M^(-1/2)+adjoint, H'=H,
%   and exp(-1i*T*H) is Euclidean unitary, while the physical
%   V=M^(-1/2)*U*M^(1/2) is M-unitary, NOT Euclidean unitary
%   for a nonconstant mass metric.
%
%   Chiral grading Gamma=diag(+I0,-I1,+I2) anticommutes with H:
%   Gamma*H*Gamma=-H. Thus H eigenvalues pair as +/-lambda,
%   U eigenphases pair as exp(-it lambda),exp(+it lambda),
%   and Gamma*U*Gamma=U'. These are STRUCTURAL invariants,
%   not a closed eigenphase law for a general triangulation.
%   The exact number of ZERO H eigenvalues equals b0+b1+b2
%   (cochain Betti numbers); there are at least this many
%   unit eigenvalues of U, with possible EXTRA resonance
%   e^{-it lambda}=1 when t lambda is a multiple of 2*pi.
%   Euler characteristic b0-b1+b2 = nV-nE+nF.
%
%   Default triangulation: oriented closed OCTAHEDRON sphere, with
%   Euler characteristic 2, (b0,b1,b2)=(1,0,1), two exact zero modes.
%
%   MODE='matrix' (default): explicit dense exact expm, N<=256.
%        'model' : dense U and physical M-unitary V, sparse H,
%                  sparse cochain d, M, D, incidence, Betti, chiral.
%        'generator': sparse H, M, D, incidence data, Betti/geometry
%                     without evaluating expm, up to N<=20000.
%        'cayley_operator': cached sparse rational unitary step
%                     [(I+itH/(2s))\(I-itH/(2s))]^s with s=STEPS
%                     (default 16), NOT exactly exp(-itH);
%                     local error O(s^-2) for fixed finite H,t.
%                     Forward, adjoint, transpose, size, multi-RHS.
%        'spectral': finite eigenphase law GENERIC UNKNOWN but
%                    symmetry and exact *guaranteed* 1-eigenvalue
%                    multiplicity lower bound b0+b1+b2.
%   Triangles may include a boundary (absolute cochains), but
%   must have consistent oriented face adjacency, each edge
%   incident to one or two faces, and nonzero triangle area.
%   All supplied vertices must participate in faces.
%
%   [176] Arnold–Falk–Winther, Acta Numerica 15 (2006), 1–155.
%         https://doi.org/10.1017/S0962492906210018
%   [177] Butt, Catterall, Pradhan & Toga, Phys. Rev. D 104 (2021), 094504.
%         https://doi.org/10.1103/PhysRevD.104.094504
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
defaultVertices=isempty(vertices);
if defaultVertices
    vertices=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
end
if nargin<2||isempty(faces)
    if ~defaultVertices
        error('simplicial_hodge_dirac:Faces', ...
            'FACES must be supplied for custom vertex positions.');
    end
    faces=[5 1 3;5 3 2;5 2 4;5 4 1; ...
           6 3 1;6 2 3;6 4 2;6 1 4];
end
if nargin<3||isempty(t),t=.7;end
if nargin<4||isempty(mode),mode='matrix';end
if nargin<5||isempty(steps),steps=16;end
if ~(isnumeric(t)&&isscalar(t)&&isreal(t)&&isfinite(t))
    error('simplicial_hodge_dirac:Time','T must be real finite scalar.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('simplicial_hodge_dirac:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','model','generator', ...
        'cayley_operator','spectral'}))
    error('simplicial_hodge_dirac:Mode', ...
        'Use matrix, model, generator, cayley_operator, spectral.');
end
if ~(isnumeric(steps)&&isscalar(steps)&&isreal(steps) ...
        &&isfinite(steps)&&steps==fix(steps)&&steps>=1&&steps<=4096)
    error('simplicial_hodge_dirac:Steps', ...
        'STEPS must be integer 1..4096.');
end
if ~(isnumeric(vertices)&&isreal(vertices)&&ismatrix(vertices) ...
        &&size(vertices,2)==3&&size(vertices,1)>=3 ...
        &&all(isfinite(vertices(:))))
    error('simplicial_hodge_dirac:Vertices', ...
        'Vertices must be a finite real nV-by-3 array.');
end
V=double(vertices);nv=size(V,1);
if ~(isnumeric(faces)&&isreal(faces)&&ismatrix(faces) ...
        &&size(faces,2)==3&&size(faces,1)>=1 ...
        &&all(isfinite(faces(:)))&&all(faces(:)>=1) ...
        &&all(faces(:)<=nv)&&all(faces(:)==fix(faces(:))))
    error('simplicial_hodge_dirac:Faces', ...
        'Faces must be integer nF-by-3 vertex indices.');
end
F=double(faces);nf=size(F,1);
if any(F(:,1)==F(:,2)||F(:,2)==F(:,3)||F(:,3)==F(:,1)) ...
        ||size(unique(sort(F,2),'rows'),1)~=nf ...
        ||numel(unique(F(:)))~=nv
    error('simplicial_hodge_dirac:Topology', ...
        'Triangles must be distinct, nondegenerate in indices, and use all vertices.');
end
a=V(F(:,1),:);b=V(F(:,2),:);c=V(F(:,3),:);
normal=cross(b-a,c-a,2);
areas=.5*sqrt(sum(normal.^2,2));
edgeScale=max(1,max(abs(V(:))));
if any(~isfinite(areas)||areas<=1e-12*edgeScale^2)
    error('simplicial_hodge_dirac:Geometry', ...
        'All faces must have positive nondegenerate area.');
end
faceEdges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
[edges,~,edgeIndex]=unique(sort(faceEdges,2),'rows');
ne=size(edges,1);N=nv+ne+nf;
if N>20000
    error('simplicial_hodge_dirac:Order', ...
        'Total 0/1/2 cochain dimension must be <=20000.');
end
incidenceCount=accumarray(edgeIndex,1,[ne 1]);
if any(incidenceCount>2)
    error('simplicial_hodge_dirac:Manifold', ...
        'An edge may belong to at most two triangles.');
end
edgeSign=2*(faceEdges(:,1)<faceEdges(:,2))-1;
D0=sparse([(1:ne)';(1:ne)'], ...
    [edges(:,1);edges(:,2)],[-ones(ne,1);ones(ne,1)],ne,nv);
D1=sparse(repmat((1:nf)',3,1),edgeIndex,edgeSign,nf,ne);
internal=find(incidenceCount==2);
if any(abs(full(sum(D1(:,internal),1)))>0)
    error('simplicial_hodge_dirac:Orientation', ...
        'Adjacent faces must induce opposite oriented shared edges.');
end
if nnz(D1*D0)>0
    error('simplicial_hodge_dirac:Internal', ...
        'Oriented incidence must satisfy d1*d0=0.');
end
edgeLengths=sqrt(sum((V(edges(:,2),:)-V(edges(:,1),:)).^2,2));
if any(edgeLengths<=1e-12*edgeScale)
    error('simplicial_hodge_dirac:Geometry','Mesh has zero-length edges.');
end
m0=accumarray(F(:),repmat(areas/3,3,1),[nv 1]);
m1=accumarray(edgeIndex, ...
    repmat(areas/3,3,1)./(edgeLengths(edgeIndex).^2),[ne 1]);
m2=1./areas;
mass=[m0;m1;m2];
if any(~isfinite(mass)||mass<=0)
    error('simplicial_hodge_dirac:Metric', ...
        'Geometric diagonal metric must be finite and positive.');
end
d=spalloc(N,N,2*nnz(D0)+2*nnz(D1));
d(nv+(1:ne),1:nv)=D0;
d(nv+ne+(1:nf),nv+(1:ne))=D1;
R=spdiags(sqrt(mass),0,N,N);
Ri=spdiags(1./sqrt(mass),0,N,N);
M=spdiags(mass,0,N,N);
A=R*d*Ri;
H=A+A';
% Combinatorial Betti numbers on a connected/disconnected orientable
% simplicial 2-manifold, with possible boundary, computed without
% numerical rank sensitivity to triangulation refinement.
components=conncomp(graph(edges(:,1),edges(:,2),[],nv));
b0=max(components);
boundaryEdge=find(incidenceCount==1);
boundaryComponents=unique(components(edges(boundaryEdge,1)));
b2=b0-numel(boundaryComponents);
chi=nv-ne+nf;
b1=b0+b2-chi;
if b1<0
    error('simplicial_hodge_dirac:Topology', ...
        'Inconsistent two-manifold Euler/Betti invariants.');
end
betti=[b0 b1 b2];
chirality=[ones(nv,1);-ones(ne,1);ones(nf,1)];
common=struct('H',H,'metric',M,'d',d,'D0',D0,'D1',D1, ...
    'edges',edges,'faces',F,'vertices',V, ...
    'face_areas',areas,'edge_lengths',edgeLengths, ...
    'cochain_dimensions',[nv ne nf], ...
    'dimension',N,'betti',betti,'euler_characteristic',chi, ...
    'chiral_grading',chirality, ...
    'boundary_edge_count',numel(boundaryEdge), ...
    'notes',['Real symmetric Euclidean Hodge–Dirac H with explicit ', ...
        'positive LUMPED cochain metric. Not spin-Dirac; ', ...
        'metric is not the exact Whitney Hodge star. ', ...
        'Zero modes equal total Betti number, chiral spectrum +/- pairs.']);
if strcmp(mode,'generator'),U=common;return;end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('simplicial_hodge_dirac',N);
    U.finite.kind='unknown';
    U.known_zero_modes=sum(betti);
    U.guaranteed_unit_eigenvalues=sum(betti);
    U.betti=betti;U.euler_characteristic=chi;
    U.chiral_pairing='eigenphases occur in conjugate pairs';
    U.notes=['Generic finite eigenspectrum UNKNOWN, but cochain ', ...
        'Betti numbers give guaranteed exp(-it*0)=1 phase modes. ', ...
        'Additional unit phases may arise from time aliasing. ', ...
        'This is Hodge–de Rham Kähler–Dirac, NOT spin-Dirac.'];
    U.reference='Arnold--Falk--Winther (2006); Catterall (2021)';
    return
end
if strcmp(mode,'cayley_operator')
    h=t/steps;
    Ip=speye(N)+.5i*h*H;
    Im=speye(N)-.5i*h*H;
    plusFactor=decomposition(Ip,'lu');
    minusFactor=decomposition(Im,'lu');
    apply=@(X,flag)um_simplicial_hodge_cayley( ...
        X,flag,N,steps,Ip,Im,plusFactor,minusFactor);
    U=struct('apply',apply,'dimension',N,'steps',steps, ...
        'time',t,'approximation', ...
        'unitary_Cayley_rational_not_exact_expm', ...
        'betti',betti,'generator',H, ...
        'notes',['Exact-in-arithmetic unitary rational Cayley step ', ...
        'with phase pairing and metric-normalization; ', ...
        'NOT the exact exp(-itH), fixed-t error O(steps^-2).']);
    return
end
if N>256
    error('simplicial_hodge_dirac:DenseOrder', ...
        'Dense expm matrix/model requires cochain dimension <=256.');
end
Q=expm(-1i*t*full(H));
if strcmp(mode,'matrix'),U=Q;return;end
% Physical cochain basis is M-weighted; the normalized matrix Q
% is unitary in the Euclidean inner product instead.
common.U=Q;
common.D_physical=Ri*H*R;
common.V_physical=Ri*Q*R;
common.time=t;
U=common;
end
