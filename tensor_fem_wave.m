function [U,properties]=tensor_fem_wave(grid,L,c,omega0,h,method,mode)
%TENSOR_FEM_WAVE Separable 2D/3D finite-element wave orthogonal steps.
%   U=TENSOR_FEM_WAVE(GRID,L,C,OMEGA0,H,METHOD,MODE), where GRID is a
%   two- or three-component vector of positive interior Dirichlet-node
%   counts, returns a real orthogonal 2N-by-2N energy-normalized wave
%   propagator, N=prod(GRID), on a rectangular 2D/3D physical domain.
%   GRID=[nx ny] or [nx ny nz], x is MATLAB's FASTEST tensor index.
%   The mesh consists of conforming bilinear/trilinear tensor products
%   of continuous piecewise-linear 1D FE functions (consistent mass).
%   L and C are positive vectors giving domain side lengths and
%   coordinate wave speeds. OMEGA0>=0 is the uniform onsite frequency.
%
%   The physical equation is Mass*q''+K*q=0 with
%       Mass=M_d kron ... kron M_1,
%       K=OMEGA0^2*Mass+sum_j C(j)^2 * (M_d kron...
%          kron K_j kron ... kron M_1).
%   M_j=(dx_j/6)*tridiag(1,4,1) and
%   K_j=(1/dx_j)*tridiag(-1,2,-1) are CONSISTENT 1D FE matrices.
%   Each 1D Dirichlet sine mode has mass m_j=dx/3*(2+cos(theta)),
%   stiffness k_j=2/dx*(1-cos(theta)), theta=pi*j/(n+1).
%   Exact multi-index frequencies:
%       w^2=OMEGA0^2+sum_j C(j)^2*k_j/m_j.
%
%   METHOD: 'midpoint' (default), 'gauss2', 'gauss3' (Gauss-Legendre
%   Padé time integrators), 'verlet' (velocity/Stoermer-Verlet).
%   Midpoint/Gauss preserve the original positive quadratic energy;
%   stable Verlet preserves a positive MODIFIED quadratic energy only
%   when H*max(w)<2; this strict inequality is validated.
%   The returned real orthogonal U is also symplectic. The physical
%   state [q;p], p=Mass*q', has a generally nonorthogonal symplectic
%   step S=T\U*T with S'*G*S=G, G=T'*T>0.
%
%   MODE='matrix' (default): dense 2N-by-2N U (N<=128).
%   MODE='operator': fast O(N log N) multidimensional orthonormal
%       DST-I action O(X,flag), no dense matrices, N<=2e6.
%   MODE='spectral': exact finite trace eigenphases and moments
%       from all 2N conjugate-paired Gauss/Verlet phase samples.
%       NO universal limiting density / van Hove law is asserted.
%   MODE='pencil': physical S, positive metric G, normalization T,
%       sparse consistent Mass and K, and A=G*S,B=G (N<=128).
%       B\A=S is circle-spectral but not generally orthogonal.
%   MODE='polynomial' (VERLET only): real symmetric palindromic
%       quadratic P(z)=z^2*Mass+z*(H^2*K-2*Mass)+Mass.
%       This has unit-circle roots only for strict stable Verlet.
%
%   References:
%   [101] Hairer, Lubich, Wanner, Geometric Numerical Integration (2006).
%   https://doi.org/10.1007/3-540-30666-8
%   [103] Cockburn, Du, Sanchez, Front. Appl. Math. Stat. 9 (2023).
%   https://doi.org/10.3389/fams.2023.1165371
%   [141] Loveland, Valseth, Lukac, Dawson, J. Comput. Sci. 64 (2022).
%   https://doi.org/10.1016/j.jocs.2022.101831
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(grid)&&isvector(grid)&&isreal(grid) ...
        &&(numel(grid)==2||numel(grid)==3) ...
        &&all(isfinite(grid))&&all(grid==floor(grid)) ...
        &&all(grid>=1)&&prod(double(grid))<=2e6)
    error('tensor_fem_wave:Grid', ...
        'GRID must have 2 or 3 positive integer sizes with prod<=2e6.');
end
grid=double(grid(:).');d=numel(grid);N=prod(grid);
if nargin<2||isempty(L),L=ones(1,d);end
if nargin<3||isempty(c),c=ones(1,d);end
if nargin<4||isempty(omega0),omega0=.8;end
if ~(isnumeric(L)&&isvector(L)&&numel(L)==d&&isreal(L) ...
        &&all(isfinite(L))&&all(L>0))
    error('tensor_fem_wave:Lengths', ...
        'L must contain one positive domain length per coordinate.');
end
if ~(isnumeric(c)&&isvector(c)&&numel(c)==d&&isreal(c) ...
        &&all(isfinite(c))&&all(c>0))
    error('tensor_fem_wave:Speeds', ...
        'C must contain one positive wave speed per coordinate.');
end
L=double(L(:).');c=double(c(:).');
validateattributes(omega0,{'numeric'},{'scalar','real','finite','nonnegative'});
if nargin<5||isempty(h),h=.12*min((L./(grid+1))./c);end
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
if nargin<6||isempty(method),method='midpoint';end
if nargin<7||isempty(mode),mode='matrix';end
if ~(ischar(method)||(isstring(method)&&isscalar(method))) ...
        || ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('tensor_fem_wave:Mode','METHOD and MODE must be scalar text.');
end
method=lower(char(method));mode=lower(char(mode));
if ~any(strcmp(method,{'midpoint','gauss2','gauss3','verlet'}))
    error('tensor_fem_wave:Method','Unknown time integration method.');
end
if ~any(strcmp(mode,{'matrix','operator','spectral','pencil','polynomial'}))
    error('tensor_fem_wave:Mode','Unknown output mode.');
end
if strcmp(mode,'polynomial')&&~strcmp(method,'verlet')
    error('tensor_fem_wave:Polynomial', ...
        'Only Verlet has the implemented palindromic wave recurrence.');
end
if any(strcmp(mode,{'matrix','pencil'}))&&N>128
    error('tensor_fem_wave:ExplicitOrder', ...
        'Dense matrix and physical pencil require prod(GRID)<=128.');
end
if strcmp(mode,'polynomial')&&N>10000
    error('tensor_fem_wave:ExplicitOrder', ...
        'Explicit sparse polynomial requires prod(GRID)<=10000.');
end
% Compute separable exact consistent-FE generalized sine symbols.
mass1=cell(d,1);stiff1=cell(d,1);ratios=cell(d,1);
for j=1:d
    dx=L(j)/(grid(j)+1);
    theta=pi*(1:grid(j))'/(grid(j)+1);
    mass1{j}=(dx/3)*(2+cos(theta));
    stiff1{j}=(2/dx)*(1-cos(theta));
    ratios{j}=c(j)^2*stiff1{j}./mass1{j};
end
massSpectral=ones(grid);
w2=omega0^2*ones(grid);
for j=1:d
    shape=ones(1,d);shape(j)=grid(j);
    massSpectral=massSpectral.*reshape(mass1{j},shape);
    w2=w2+reshape(ratios{j},shape);
end
massSpectral=massSpectral(:);w=sqrt(w2(:));
stiffSpectral=massSpectral.*w2(:);
if strcmp(method,'verlet')
    if h*max(w)>=2
        error('tensor_fem_wave:Stability', ...
          'Positive Verlet metric requires strict h*max(w)<2.');
    end
    x=h*w;
    z=(1-.5*x.^2)-1i*x.*sqrt(1-.25*x.^2);
else
    switch method
        case 'midpoint',coeff=[1 .5];
        case 'gauss2',  coeff=[1 .5 1/12];
        case 'gauss3',  coeff=[1 .5 1/10 1/120];
    end
    z=polyval(fliplr(coeff),-1i*h*w)./ ...
      polyval(fliplr(coeff),1i*h*w);
end
if strcmp(mode,'spectral')
    U=um_geometric_spectral('tensor_fem_wave',z,w, ...
         massSpectral,stiffSpectral);
    U.grid_shape=grid;
    U.spatial_dimension=d;
    U.domain_lengths=L;
    U.wave_speeds=c;
    U.method=method;
    U.time_step=h;
    U.notes=['Exact finite tensor-product consistent-FE generalized ', ...
       'frequencies and conjugate eigenphases; repeated/degenerate ', ...
       'phases retain their multiplicities in the trace measure. ', ...
       'No analytic N->infinity or van Hove-type law is asserted.'];
    return
end
if strcmp(mode,'operator')
    U=@(X,flag)um_tensor_fem_apply(X,flag,grid,z);
    return
end
% Actual Galerkin tensors. Their Kronecker order agrees with
% MATLAB column-major physical index 1 + (ix-1) + nx*(iy-1) + ...
M1=cell(d,1);K1=cell(d,1);
for j=1:d
    dx=L(j)/(grid(j)+1);
    a=ones(grid(j),1);
    M1{j}=(dx/6)*spdiags([a 4*a a],-1:1,grid(j),grid(j));
    K1{j}=(1/dx)*spdiags([-a 2*a -a],-1:1,grid(j),grid(j));
end
Mass=1;
for j=d:-1:1,Mass=kron(Mass,M1{j});end
Stiffness=omega0^2*Mass;
for j=1:d
    term=1;
    for k=d:-1:1
        if k==j,block=K1{k};else,block=M1{k};end
        term=kron(term,block);
    end
    Stiffness=Stiffness+c(j)^2*term;
end
if strcmp(mode,'polynomial')
    U=struct('A0',Mass,'A1',h*h*Stiffness-2*Mass, ...
       'A2',Mass,'mass',Mass,'stiffness',Stiffness, ...
       'dimension',N,'grid_shape',grid,'method',method, ...
       'classification','real_symmetric_palindromic_qep');
    return
end
% Dense exact modal representation, V=V_d kron ... kron V_1.
V=1;
for j=d:-1:1
    rows=(1:grid(j))';
    Vi=sqrt(2/(grid(j)+1))*sin(pi*(rows*rows')/(grid(j)+1));
    V=kron(V,Vi);
end
Z=(V.*z.')*V';
C=real(Z);S=-imag(Z);
U=[C,S;-S,C];
if strcmp(mode,'matrix'),return;end
if strcmp(method,'verlet')
    metricQ=stiffSpectral.*(1-.25*h*h*w.^2);
else
    metricQ=stiffSpectral;
end
metricP=1./massSpectral;
Tq=(V.*sqrt(metricQ).')*V';
Tp=(V.*sqrt(metricP).')*V';
T=blkdiag(Tq,Tp);
G=blkdiag((V.*metricQ.')*V',(V.*metricP.')*V');
Mstep=T\(U*T);
U=struct('A',G*Mstep,'B',G,'M',Mstep, ...
    'G',G,'T',T,'U',U,'mass',Mass,'stiffness',Stiffness, ...
    'grid_shape',grid,'method',method, ...
    'frequencies',w,'classification','circle_spectrum_only');
end
