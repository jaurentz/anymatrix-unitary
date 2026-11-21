function [U,properties]=flat_torus_spin_dirac(nx,ny,t,mass,twist,weights,mode)
%FLAT_TORUS_SPIN_DIRAC Exact unitary centered lattice spin-Dirac propagator.
%   U=FLAT_TORUS_SPIN_DIRAC(NX,NY,T,MASS,TWIST,WEIGHTS,MODE) defines
%   a periodic NX-by-NY square-cell discretization of the flat 2-torus,
%   with two complex spinor components per lattice site. Its real
%   dimensions are NX,NY>=3, N=NX*NY, matrix order 2*N. This is a
%   finite-difference NAIVE DIRAC lattice model: not the triangle-mesh
%   intrinsic/extrinsic Dirac discretization of Ye et al. (2018).
%
%   T is the evolution time (default .7), MASS the Dirac mass
%   (default .2). TWIST=[beta_x,beta_y] specifies spinorial boundary
%   holonomy in radians (default [0 0]): each beta=0 or pi is a
%   periodic/antiperiodic torus spin structure; arbitrary angles
%   represent additional flat U(1) twists, NOT new spin structures.
%   Lattice spacings are fixed hx=hy=1 (no claims about continuum
%   rescaling). Site index is j=1+x+NX*y, x=0,...,NX-1, y=0,...,NY-1,
%   with up/down spinor coordinates (2*j-1,2*j).
%
%   The unitary covariant shifts act as (Tx f)(x,y)=f(x+1,y)
%   and (Ty f)(x,y)=f(x,y+1), with exp(i beta_x),exp(i beta_y)
%   on their respective positive-wrap links. Put
%     Px=(Tx-Tx')/(2i), Py=(Ty-Ty')/(2i),  H=Px*sigma_x+Py*sigma_y
%         +MASS*sigma_z (tensor factors in site/spin ordering).
%   H is Hermitian and sparse (O(N) entries); U=exp(-i*T*H)
%   is exactly unitary as a finite mathematical matrix.
%   The MODE='operator' action diagonalizes H by TWISTED FFTs:
%     hx(k)=sin((2*pi*k+beta_x)/NX),
%     hy(j)=sin((2*pi*j+beta_y)/NY),
%     E(k,j)=sqrt(hx(k)^2+hy(j)^2+MASS^2).
%   Every momentum sector has eigenvalues +/-E; hence the exact
%   propagator eigenphase multiset is exp(-i*T*(+/-E)).
%   WARNING: naive centered differences have species/fermion doubling
%   at Brillouin-zone points; this is not a Wilson/staggered correction.
%
%   WEIGHTS is an OPTIONAL positive scalar or N-vector of lumped
%   site metric weights. Define M=diag(repelem(WEIGHTS,2));
%   D_phys=M^(-1/2)*H*M^(1/2), so D_phys'*M=M*D_phys, and
%   V_phys=M^(-1/2)*U*M^(1/2) is M-unitary, not necessarily
%   Euclidean unitary. The public U always uses Euclidean
%   orthonormal coordinates. WEIGHTS are an explicit coordinate
%   metric model, NOT assembled from arbitrary triangle-mesh FEM.
%
%   MODE='matrix' (default): explicit dense U for N<=196.
%        'operator': O(N log N) FFT forward/adjoint/transpose
%                    application to multiple RHS, and 'size'.
%        'model': explicit U, sparse Hermitian H, sparse metric M,
%                 sparse D_phys, M-unitary physical V_phys (N<=196).
%        'spectral': exact finite trace eigenphases/moments
%                    from sine dispersion, no dense matrix.
%   Mode operator can handle N<=1e6 and does NOT require dense
%   matrix exponentials. Generic arbitrary mesh edges, curved
%   spin connections, nonflat gauge fields and zero-mode index
%   theory are beyond this baseline.
%
%   [169] A. Donis Vela and C. W. J. Beenakker (2025),
%   https://doi.org/10.1088/1751-8121/adbac7
%   The cited work uses a DIFFERENT tangent fermion discretization,
%   not the naive centered formula implemented here.
%   [170] Ye et al. (2018), mesh Dirac operators: CONTEXT ONLY.
%   https://doi.org/10.1111/cgf.13494
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(nx)&&isscalar(nx)&&isreal(nx)&&isfinite(nx) ...
    &&nx==fix(nx)&&nx>=3&&nx<=1000)
    error('flat_torus_spin_dirac:Order','NX must be integer 3..1000.');
end
if nargin<2||isempty(ny),ny=nx;end
if ~(isnumeric(ny)&&isscalar(ny)&&isreal(ny)&&isfinite(ny) ...
    &&ny==fix(ny)&&ny>=3&&ny<=1000)
    error('flat_torus_spin_dirac:Order','NY must be integer 3..1000.');
end
n=nx*ny;
if nargin<3||isempty(t),t=.7;end
if nargin<4||isempty(mass),mass=.2;end
if nargin<5||isempty(twist),twist=[0;0];end
if nargin<6||isempty(weights),weights=1;end
if nargin<7||isempty(mode),mode='matrix';end
if ~(isnumeric(t)&&isscalar(t)&&isreal(t)&&isfinite(t) ...
    &&isnumeric(mass)&&isscalar(mass)&&isreal(mass)&&isfinite(mass))
    error('flat_torus_spin_dirac:Parameters', ...
        'T and MASS must be real finite scalars.');
end
if ~(isnumeric(twist)&&isreal(twist)&&isvector(twist) ...
    &&numel(twist)==2&&all(isfinite(twist(:))))
    error('flat_torus_spin_dirac:Twist', ...
        'TWIST must be a two-component finite real vector in radians.');
end
twist=double(twist(:));
twist=mod(twist+pi,2*pi)-pi;
if ~(isnumeric(weights)&&isreal(weights)&&isvector(weights) ...
    &&all(isfinite(weights(:)))&&all(weights(:)>0) ...
    &&(isscalar(weights)||numel(weights)==n) ...
    &&max(weights(:))<=1e12 && min(weights(:))>=1e-12)
    error('flat_torus_spin_dirac:Weights', ...
        'WEIGHTS must be positive scalar or NX*NY vector, in [1e-12,1e12].');
end
if isscalar(weights),weights=repmat(double(weights),n,1);
else,weights=double(weights(:));end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('flat_torus_spin_dirac:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','model','spectral'}))
    error('flat_torus_spin_dirac:Mode', ...
        'MODE must be matrix, operator, model, spectral.');
end
if any(strcmp(mode,{'matrix','model'}))&&n>196
    error('flat_torus_spin_dirac:ExplicitOrder', ...
        'Dense explicit/model output requires NX*NY<=196.');
end
[kx,ky]=ndgrid(0:nx-1,0:ny-1);
px=sin((2*pi*kx+twist(1))/nx);
py=sin((2*pi*ky+twist(2))/ny);
E=sqrt(px.^2+py.^2+mass^2);
if strcmp(mode,'spectral')
    lam=[exp(-1i*t*E(:));exp(1i*t*E(:))];
    U=unitary_spectral_reference('flat_torus_spin_dirac',2*n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()lam;
    U.finite.moments=@(k)torus_dirac_moments(lam,k);
    U.reference=['A. Donis Vela and C. W. J. Beenakker (2025), ', ...
        'doi:10.1088/1751-8121/adbac7 (different lattice scheme)'];
    U.notes=['Exact finite naive centered-lattice sine dispersion; ', ...
        'periodic/antiperiodic torus spin structures (beta=0/pi); ', ...
        'other twists represent flat U(1) holonomy. Fermion ', ...
        'doubling is present; continuum Dirac eigenvalues differ. ', ...
        'Site-metric weights give M-selfadjoint physical D but ', ...
        'do not alter the Euclidean-normalized spectrum.'];
    U.lattice_size=[nx ny];U.mass=mass;U.time=t;U.twist=twist;
    return
end
apply=@(X,flag)um_flat_torus_spin_dirac_apply( ...
    X,flag,nx,ny,t,mass,twist,px,py,E);
if strcmp(mode,'operator')
    U=apply;
    return
end
U=apply(eye(2*n),'notransp');
if strcmp(mode,'matrix'),return;end

% Distinguish finite Euclidean-unitary U from the positive-metric
% physical propagation V=M^(-1/2)*U*M^(1/2).
tx=sparse(1:nx,[2:nx 1], ...
    [ones(1,nx-1),exp(1i*twist(1))],nx,nx);
ty=sparse(1:ny,[2:ny 1], ...
    [ones(1,ny-1),exp(1i*twist(2))],ny,ny);
Ax=kron(speye(ny),(tx-tx')/(2i));
Ay=kron((ty-ty')/(2i),speye(nx));
sx=sparse([0,1;1,0]);
sy=sparse([0,-1i;1i,0]);
sz=sparse([1,0;0,-1]);
H=kron(Ax,sx)+kron(Ay,sy)+mass*kron(speye(n),sz);
w=repelem(weights,2);
M=spdiags(w,0,2*n,2*n);
R=spdiags(sqrt(w),0,2*n,2*n);
Ri=spdiags(1./sqrt(w),0,2*n,2*n);
D=Ri*H*R;
V=Ri*U*R;
U=struct('U',U,'H',H,'metric',M,'D_physical',D, ...
    'V_physical',V,'lattice_size',[nx ny], ...
    'mass',mass,'time',t,'twist',twist, ...
    'site_weights',weights,'momentum_energy',E, ...
    'notes',['Normalized Euclidean unitary U=exp(-i*t*H). ', ...
      'The physical D=M^(-1/2)*H*M^(1/2) is M-selfadjoint, ', ...
      'V=M^(-1/2)*U*M^(1/2) is M-unitary, and these ', ...
      'arbitrary positive site weights are not a curved ', ...
      'triangle-mesh FEM mass. Naive finite-difference doubling.']);
end

function m=torus_dirac_moments(z,k)
validateattributes(k,{'numeric'},{'integer','real','finite'});
m=reshape(arrayfun(@(p)mean(z.^p),k(:)),size(k));
end
