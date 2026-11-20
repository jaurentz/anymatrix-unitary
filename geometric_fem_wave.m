function [U,properties] = geometric_fem_wave(n,L,c,omega0,h,method,mode)
%GEOMETRIC_FEM_WAVE Unitary transforms of a 1D finite-element wave step.
%   U=GEOMETRIC_FEM_WAVE(N,L,C,OMEGA0,H,METHOD) discretizes
%       q_tt-c^2 q_xx+OMEGA0^2 q=0 on (0,L), zero Dirichlet endpoints,
%   by uniform continuous piecewise-linear finite elements with N interior
%   degrees of freedom, and advances the canonical Hamiltonian system
%       q'=Mass\p, p'=-K*q
%   using implicit midpoint, Gauss-Legendre order four/six or velocity
%   Verlet. The returned 2N-by-2N matrix is orthogonal-symplectic.
%   Default L=1, C=1, OMEGA0=1, METHOD='midpoint', H=.1*L/(N+1)/C.
%   METHODS: 'midpoint', 'gauss2', 'gauss3', 'verlet'.
%
%   'pencil' returns struct with A=G*Mstep, B=G, physical Mstep=M,
%   metric G>0, normalization T, mass/stiffness, and normalized U.
%   This pencil has unit-circle eigenvalues, but B\A is not generally
%   Euclidean unitary. For Verlet, 'polynomial' returns the symmetric
%   palindromic quadratic eigenvalue problem
%       z^2*Mass+z*(H^2*K-2*Mass)+Mass.
%   'spectral' returns the exact 2N eigenvalues from generalized
%   finite-element sine-mode frequencies, without eig.
%   'operator' returns O(X,FLAG), a fast DST-I matrix-free action.
%
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration,
%   2nd ed. (2006). [97] https://doi.org/10.1007/3-540-30666-8
%   B. Cockburn, S. Du, M. A. Sanchez, Front. Appl. Math. Stat. 9 (2023).
%   [99] https://doi.org/10.3389/fams.2023.1165371
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(L),L=1;end
if nargin<3 || isempty(c),c=1;end
if nargin<4 || isempty(omega0),omega0=1;end
if nargin<5 || isempty(h),h=0.1*L/(n+1)/c;end
if nargin<6 || isempty(method),method='midpoint';end
if nargin<7 || isempty(mode),mode='matrix';end
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(omega0,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
if ~(ischar(method)||isstring(method)) || ~(ischar(mode)||isstring(mode))
    error('geometric_fem_wave:Mode','METHOD and MODE must be text.');
end
method=lower(char(method));mode=lower(char(mode));
if ~any(strcmp(method,{'midpoint','gauss2','gauss3','verlet'}))
    error('geometric_fem_wave:Method','Unknown integration method.');
end
if ~any(strcmp(mode,{'matrix','operator','pencil','spectral','polynomial'}))
    error('geometric_fem_wave:Mode','Unknown output mode.');
end
if strcmp(mode,'polynomial') && ~strcmp(method,'verlet')
    error('geometric_fem_wave:Polynomial', ...
        'Position-recurrence polynomial mode is implemented for Verlet.');
end
dx=L/(n+1);
theta=pi*(1:n)'/(n+1);
massvals=(dx/3)*(2+cos(theta));
stiffvals=(2*c*c/dx)*(1-cos(theta))+omega0^2*massvals;
w=sqrt(stiffvals./massvals);
if strcmp(method,'verlet')
    if h*max(w)>=2
        error('geometric_fem_wave:Stability', ...
            'Verlet requires h*max(finite-element frequencies)<2.');
    end
    x=h*w;z=1-.5*x.^2-1i*x.*sqrt(1-.25*x.^2);
else
    switch method
        case 'midpoint',coeff=[1 .5];
        case 'gauss2',coeff=[1 .5 1/12];
        case 'gauss3',coeff=[1 .5 1/10 1/120];
    end
    z=polyval(fliplr(coeff),-1i*h*w)./ ...
      polyval(fliplr(coeff),1i*h*w);
end
if strcmp(mode,'spectral')
    U=struct('dimension',2*n,'kind','exact', ...
        'eigenvalues',[z;conj(z)],'frequencies',w, ...
        'method',method,'mass_eigenvalues',massvals, ...
        'stiffness_eigenvalues',stiffvals);
    return
end
if strcmp(mode,'operator')
    U=@(X,flag)um_geometric_wave_apply(X,flag,z);
    return
end
a=ones(n,1);
Mass=(dx/6)*spdiags([a 4*a a],-1:1,n,n);
Lapl=(1/dx)*spdiags([-a 2*a -a],-1:1,n,n);
K=c*c*Lapl+omega0^2*Mass;
if strcmp(mode,'polynomial')
    U=struct('A0',Mass,'A1',h*h*K-2*Mass,'A2',Mass, ...
        'mass',Mass,'stiffness',K,'dimension',n,'method',method);
    return
end
j=(1:n)';
V=sqrt(2/(n+1))*sin(pi*(j*j')/(n+1));
W=(V.*z.')*V';
C=real(W);S=-imag(W);
U=[C,S;-S,C];
if strcmp(mode,'matrix'),return;end
if strcmp(method,'verlet')
    metricQ=stiffvals.*(1-.25*h*h*w.^2);
else
    metricQ=stiffvals;
end
metricP=1./massvals;
Tq=(V.*sqrt(metricQ).')*V';
Tp=(V.*sqrt(metricP).')*V';
T=blkdiag(Tq,Tp);
G=blkdiag((V.*metricQ.')*V',(V.*metricP.')*V');
Mstep=T\(U*T);
A=G*Mstep;B=G;
U=struct('A',A,'B',B,'M',Mstep,'G',G,'T',T, ...
    'U',U,'mass',Mass,'stiffness',K,'method',method, ...
    'classification','circle_spectrum_only');
end
