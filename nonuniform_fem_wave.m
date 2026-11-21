function [U,properties]=nonuniform_fem_wave(n,L,c,omega0,h,method,mode,opts)
%NONUNIFORM_FEM_WAVE Noncommuting 1D consistent-FE wave Gauss propagators.
%   U=NONUNIFORM_FEM_WAVE(N,L,C,OMEGA0,H,METHOD,MODE,OPTS)
%   returns a real orthogonal 2N-by-2N energy-normalized midpoint /
%   Gauss-Legendre wave step for genuinely nonuniform or heterogeneous
%   1D continuous piecewise-linear finite elements, homogeneous Dirichlet
%   boundaries. The physical state x=[q;p], p=Mass*q', obeys
%      Mass*q''+K*q=0, A0=[0,Mass\I;-K,0], x'=A0*x,
%      G=diag(K,Mass\I), A0'*G+G*A0=0, A0'*Jc+Jc*A0=0.
%   Mass and K are sparse real SPD; generically Mass*K ~= K*Mass.
%   The physical Gauss step S=R_s(h*A0) is symplectic and G-orthogonal.
%   T'*T=G, U=T*S/T is real orthogonal but GENERALLY NOT standard-Jc
%   symplectic. Its symplectic form is Jtilde=T'\Jc/T.
%
%   N>=2, L>0 domain length, C>0 nominal wave speed, OMEGA0>0
%   onsite frequency (allows 0 with positive Dirichlet stiffness).
%   Default L=C=1, OMEGA0=.8, H=.12*L/(N+1)/C,
%   METHOD='midpoint', MODE='matrix'.
%
%   OPTS is scalar struct with optional fields:
%     mesh_amplitude  default .27, 0<=value<.9;
%     speed_amplitude default .34, 0<=value<.9;
%     element_lengths default [], else positive N+1-vector summing L;
%     element_speeds  default [], else positive N+1-vector.
%   Default element lengths proportional to
%      1+mesh_amplitude*sin(sqrt(2)*e+.4), e=1:N+1.
%   Default speeds are C*(1+speed_amplitude*cos(.81*e+.2)).
%   Both amplitudes set to 0 -> uniform constant-coefficient 1D FE
%   benchmark from GEOMETRIC_FEM_WAVE; its known exact finite phases
%   are then recovered in SPECTRAL mode. Explicit element vectors
%   override corresponding amplitudes.
%
%   MODES:
%     'matrix'      dense 2N-by-2N orthogonal step; N<=128.
%     'operator'    cached sparse generalized Cayley-factor LU solves
%                   with sparse Cholesky energy whitening, N<=5000.
%     'pencil'      A=G*S, B=G, S physical canonical symplectic step,
%                   G,T, normalized U and transformed Jtilde; N<=128.
%     'coefficients' sparse Mass,K, nodal/element geometry, commutator.
%     'quadratic'  continuous wave QEP lambda^2*Mass+K; not unitary.
%     'spectral'   finite trace phase law UNKNOWN for heterogeneous /
%                   nonuniform systems, EXACT for uniform limit only.
%
%   Each s-stage [s/s] Pade factor R_s(z)=prod_j (1+a_j*z)/(1-a_j*z)
%   is applied to A0 via sparse generalized first-order matrices
%      E=diag(Mass,I), L0=[0,I;-K,0], E*A0=L0,
%      (I-a*h*A0)\(I+a*h*A0)=(E-a*h*L0)\(E+a*h*L0).
%   No dense Mass inverse and NO invalid DST diagonalization for a
%   noncommuting heterogeneous mesh. Factor inverses supply adjoints.
%   Positive energy normalization T can be ill conditioned; an
%   ill-conditioned Mass or K is rejected, not silently whitened.
%
%   [101] E. Hairer, C. Lubich, G. Wanner (2006).
%   https://doi.org/10.1007/3-540-30666-8
%   [103] B. Cockburn, S. Du, M. A. Sanchez (2023).
%   https://doi.org/10.3389/fams.2023.1165371
properties={'square','unitary','orthogonal','normal','real', ...
            'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
      &&n==floor(n)&&n>=2&&n<=5000)
    error('nonuniform_fem_wave:Order','N must be an integer from 2 to 5000.');
end
if nargin<2||isempty(L),L=1;end
if nargin<3||isempty(c),c=1;end
if nargin<4||isempty(omega0),omega0=.8;end
if nargin<5||isempty(h),h=.12*L/(n+1)/c;end
if nargin<6||isempty(method),method='midpoint';end
if nargin<7||isempty(mode),mode='matrix';end
if nargin<8||isempty(opts),opts=struct();end
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(omega0,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
if ~(ischar(method)||(isstring(method)&&isscalar(method))) ...
  || ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('nonuniform_fem_wave:Mode','METHOD and MODE must be scalar text.');
end
method=lower(char(method));mode=lower(char(mode));
if ~any(strcmp(method,{'midpoint','gauss2','gauss3'}))
    error('nonuniform_fem_wave:Method','Unknown Gauss time integrator.');
end
if ~any(strcmp(mode,{'matrix','operator','pencil', ...
                     'coefficients','quadratic','spectral'}))
    error('nonuniform_fem_wave:Mode','Unknown output mode.');
end
if ~isstruct(opts)||~isscalar(opts)
    error('nonuniform_fem_wave:Options','OPTS must be a scalar struct.');
end
defaults=struct('mesh_amplitude',.27,'speed_amplitude',.34, ...
    'element_lengths',[],'element_speeds',[]);
names=fieldnames(opts);
for j=1:numel(names)
    if ~isfield(defaults,names{j})
        error('nonuniform_fem_wave:Options','Unrecognized mesh option.');
    end
    defaults.(names{j})=opts.(names{j});
end
p=defaults;
for name={'mesh_amplitude','speed_amplitude'}
    value=p.(name{1});
    if ~(isnumeric(value)&&isscalar(value)&&isreal(value) ...
            &&isfinite(value)&&value>=0&&value<.9)
        error('nonuniform_fem_wave:Options', ...
            'Mesh and speed amplitudes must belong to [0,.9).');
    end
end
if isempty(p.element_lengths)
    weight=1+p.mesh_amplitude*sin(sqrt(2)*(1:n+1)'+.4);
    widths=L*weight/sum(weight);
else
    widths=p.element_lengths(:);
    if ~(isnumeric(widths)&&numel(widths)==n+1 ...
            &&isreal(widths)&&all(isfinite(widths))&&all(widths>0) ...
            &&abs(sum(widths)-L)<=1e-10*L)
        error('nonuniform_fem_wave:Mesh', ...
            'Element lengths must be positive N+1 vector summing L.');
    end
end
if isempty(p.element_speeds)
    speeds=c*(1+p.speed_amplitude*cos(.81*(1:n+1)'+.2));
else
    speeds=p.element_speeds(:);
    if ~(isnumeric(speeds)&&numel(speeds)==n+1 ...
            &&isreal(speeds)&&all(isfinite(speeds))&&all(speeds>0))
        error('nonuniform_fem_wave:Mesh', ...
            'Element speeds must be a positive N+1 vector.');
    end
end
% Assemble actual continuous piecewise-linear FE element contributions.
% Global nodes 0,..,N+1; fixed Dirichlet endpoints 0 and N+1.
md=zeros(n,1);kd=zeros(n,1);
mo=zeros(n-1,1);ko=zeros(n-1,1);
for e=1:n+1
    dx=widths(e);ce=speeds(e);
    mloc=dx/3; kloc=ce^2/dx;
    a=e-1; b=e;
    if a>=1
        md(a)=md(a)+mloc;kd(a)=kd(a)+kloc;
    end
    if b<=n
        md(b)=md(b)+mloc;kd(b)=kd(b)+kloc;
    end
    if a>=1 && b<=n
        mo(a)=dx/6;ko(a)=-kloc;
    end
end
Mass=spdiags(md,0,n,n)+sparse(1:n-1,2:n,mo,n,n) ...
    +sparse(2:n,1:n-1,mo,n,n);
K=spdiags(kd,0,n,n)+sparse(1:n-1,2:n,ko,n,n) ...
    +sparse(2:n,1:n-1,ko,n,n)+omega0^2*Mass;
% Reject numerically unresolvable positive-metric normalizations.
[Rm,pm]=chol(Mass);
[Rk,pk]=chol(K);
if pm~=0||pk~=0||condest(Mass)>1e12||condest(K)>1e12
    error('nonuniform_fem_wave:Conditioning', ...
        'Mass/stiffness are not well-conditioned SPD matrices.');
end
comm=norm(Mass*K-K*Mass,'fro');
nodes=[0;cumsum(widths)];
if strcmp(mode,'coefficients')||strcmp(mode,'quadratic')
    U=struct('mass',Mass,'stiffness',K, ...
      'element_lengths',widths,'element_speeds',speeds, ...
      'nodes',nodes,'commutator_fro',comm,'dimension',n);
    if strcmp(mode,'quadratic')
        U.A0=K;U.A1=sparse(n,n);U.A2=Mass;
        U.classification='continuous_symmetric_wave_quadratic';
    end
    return
end
uniform=(max(abs(widths-L/(n+1)))<1e-12*L ...
   &&max(abs(speeds-c))<1e-12*max(1,c));
if strcmp(mode,'spectral')
    if uniform
        theta=pi*(1:n)'/(n+1);
        m=(L/(n+1)/3)*(2+cos(theta));
        k=(2*c^2/(L/(n+1)))*(1-cos(theta))+omega0^2*m;
        w=sqrt(k./m);
        coeff=um_nonuniform_wave_pade(method);
        z=polyval(fliplr(coeff),-1i*h*w)./ ...
          polyval(fliplr(coeff),1i*h*w);
        U=um_geometric_spectral('nonuniform_fem_wave',z,w,m,k);
        U.notes=['Uniform 1D limit only: exact generalized FE sine ', ...
          'frequencies and Gauss eigenphases. Heterogeneous models ', ...
          'have no asserted analytic finite trace phase law.'];
    else
        blank=struct('kind','unknown','density',[],'support',[], ...
            'atoms',[],'moments',[]);
        U=struct('family','nonuniform_fem_wave','dimension',2*n, ...
          'normalization','probability (density with respect to dtheta)', ...
          'finite',struct('kind','unknown','eigenvalues',[], ...
                          'moments',[],'atoms',[]), ...
          'limit',blank,'ensemble',blank,'underlying',blank, ...
          'reference','Hairer, Lubich and Wanner (2006); Cockburn, Du and Sanchez (2023)', ...
          'notes',['For noncommuting nonuniform or heterogeneous SPD ', ...
             'finite elements there is no claimed closed-form finite ', ...
             'trace eigenphase law; eig(U) is a numerical benchmark, ', ...
             'not an analytic reference.']);
    end
    U.mass_commutator_fro=comm;
    U.method=method;
    return
end
if any(strcmp(mode,{'matrix','pencil'}))&&n>128
    error('nonuniform_fem_wave:ExplicitOrder', ...
        'Explicit matrix/pencil requires N<=128; use operator.');
end
% Sparse generalized first-order system E*x'=L0*x, A0=E\L0.
I=speye(n);
E=blkdiag(Mass,I);
L0=[sparse(n,n),I;-K,sparse(n,n)];
coeff=um_nonuniform_wave_pade(method);
alpha=-1./roots(fliplr(coeff));
if strcmp(mode,'operator')
    direct=cell(numel(alpha),1);reverse=cell(numel(alpha),1);
    for j=1:numel(alpha)
        direct{j}=decomposition(E-alpha(j)*h*L0,'lu');
        reverse{j}=decomposition(E+alpha(j)*h*L0,'lu');
    end
    U=@(X,flag)um_nonuniform_fem_apply(X,flag,E,L0,alpha,h,Rk,Rm,direct,reverse);
    return
end
S=speye(2*n);
for j=1:numel(alpha)
    S=(E-alpha(j)*h*L0)\((E+alpha(j)*h*L0)*S);
end
S=real(full(S));
T=blkdiag(full(Rk),full(Rm'\speye(n)));
U=real(T*S/T);
if strcmp(mode,'matrix'),return;end
G=blkdiag(full(K),full(Mass\I));
Jc=[zeros(n),eye(n);-eye(n),zeros(n)];
U=struct('A',G*S,'B',G,'Mstep',S,'G',G,'T',T,'U',U, ...
    'mass',Mass,'stiffness',K,'physicalGenerator',full(E\L0), ...
    'canonicalSymplectic',Jc, ...
    'normalizedSymplectic',((T')\Jc)/T, ...
    'commutator_fro',comm,'nodes',nodes, ...
    'element_lengths',widths,'element_speeds',speeds, ...
    'method',method,'classification','circle_spectrum_only');
end
