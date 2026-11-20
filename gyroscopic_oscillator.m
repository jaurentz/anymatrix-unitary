function [U,properties]=gyroscopic_oscillator(n,h,omega0,kappa,gamma,method,mode)
%GYROSCOPIC_OSCILLATOR Sparse gyroscopic vibration QEP and unit-circle pencil.
%   U=GYROSCOPIC_OSCILLATOR(N,H,OMEGA0,KAPPA,GAMMA,METHOD,MODE)
%   gives an orthogonal 2N-by-2N energy-normalized Gauss-Legendre step
%   for M*q''+C*q'+K*q=0, with M,K real SPD and C real skew.
%   N>=2. Defaults H=.25, OMEGA0=1, KAPPA=1, GAMMA=.35.
%   METHOD: 'midpoint' (default), 'gauss2', 'gauss3'.
%   MODE: 'matrix' (default), 'operator', 'pencil', 'quadratic', 'spectral'.
%   'operator' returns O(X,FLAG) with sparse linear solves.
%   'pencil' returns physical A=G*S,B=G,Mstep=S,G,T,U, M,K,C,
%     physicalGenerator=A0 and canonicalMap=Bc.
%   'quadratic' returns A0=K,A1=C,A2=M for the continuous QEP.
%   'spectral' reports unknown generic analytic finite phase law.
%   The pencil is positive-metric circle-spectral, not generally a
%   strict Euclidean-unitary pencil.
%
%   For state x=[q;v], A0=[0,I;-M\K,-M\C], G=diag(K,M),
%   A0'*G+G*A0=0, S=R_s(h*A0), S'*G*S=G, U=T*S/T,
%   where T=chol(G). Canonical momentum p=M*v+C*q/2 and
%   canonicalMap=[I,0;C/2,M]. The normalized U does not generally
%   preserve canonical J even though the correctly transformed S does.
%
%   The heterogeneous nearest-neighbor benchmark coefficients are a
%   reproducible specialization, not measured experimental input.
%   J. Qian and W.-W. Lin, J. Sound Vib. 306 (2007), 284-296 [125].
%     https://doi.org/10.1016/j.jsv.2007.05.009
%   T. Ohsawa, SIAM J. Numer. Anal. 61 (2023), 1293-1315 [126].
%     https://doi.org/10.1137/22M1517718
%   E. Hairer, C. Lubich and G. Wanner (2006) [101].
%     https://doi.org/10.1007/3-540-30666-8
%   I. Kuzmanovic Ivicic and S. Miodragovic (2024) [127].
%     https://doi.org/10.1016/j.cam.2024.116061
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',2});
if nargin<2||isempty(h),h=.25;end
if nargin<3||isempty(omega0),omega0=1;end
if nargin<4||isempty(kappa),kappa=1;end
if nargin<5||isempty(gamma),gamma=.35;end
if nargin<6||isempty(method),method='midpoint';end
if nargin<7||isempty(mode),mode='matrix';end
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(omega0,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(kappa,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite','nonnegative'});
if ~(ischar(method)||(isstring(method)&&isscalar(method)))|| ...
        ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('gyroscopic_oscillator:Mode','METHOD and MODE must be text.');
end
method=lower(char(method)); mode=lower(char(mode));
if ~any(strcmp(method,{'midpoint','gauss2','gauss3'}))
    error('gyroscopic_oscillator:Method','Unknown method.');
end
if ~any(strcmp(mode,{'matrix','operator','pencil','quadratic','spectral'}))
    error('gyroscopic_oscillator:Mode','Unknown output mode.');
end
j=(1:n)';
Mass=spdiags(1+.2*j/n+.05*sin(j),0,n,n);
D=sparse([(1:n-1)';(1:n-1)'],[(1:n-1)';(2:n)'], ...
    [-ones(n-1,1);ones(n-1,1)],n-1,n);
K=omega0^2*Mass+kappa*(D'*D)+ ...
    .15*spdiags(1+.1*cos(j),0,n,n);
E=sparse(2:n,1:n-1,1,n,n);
C=gamma*(E-E');
if strcmp(mode,'quadratic')
    U=struct('A0',K,'A1',C,'A2',Mass,'mass',Mass, ...
        'stiffness',K,'gyroscopic',C, ...
        'classification','gyroscopic_continuous_qep');
    return
end
if strcmp(mode,'spectral')
    U=struct('family','gyroscopic_oscillator','dimension',2*n, ...
        'method',method,'finite',struct('kind','unknown', ...
        'note','No generic analytic finite phase law.'), ...
        'reference_measure','equal-weight trace measure');
    return
end
I=speye(n);
A0=[sparse(n,n),I;-(Mass\K),-(Mass\C)];
G=blkdiag(K,Mass);
T=chol(G);
Z=h*A0;I2=speye(2*n);
switch method
    case 'midpoint'
        P=I2+.5*Z;Q=I2-.5*Z;
    case 'gauss2'
        Z2=Z*Z;
        P=I2+.5*Z+Z2/12;Q=I2-.5*Z+Z2/12;
    case 'gauss3'
        Z2=Z*Z;Z3=Z2*Z;
        P=I2+.5*Z+Z2/10+Z3/120;
        Q=I2-.5*Z+Z2/10-Z3/120;
end
if strcmp(mode,'operator')
    U=@(X,flag)um_gyroscopic_oscillator_apply(X,flag,P,Q,T);
    return
end
S=Q\P;
U=full(T*S/T);
if strcmp(mode,'matrix'),return;end
canonicalMap=[I,sparse(n,n);C/2,Mass];
U=struct('A',G*S,'B',G,'Mstep',S,'G',G,'T',T,'U',U, ...
    'mass',Mass,'stiffness',K,'gyroscopic',C, ...
    'physicalGenerator',A0,'canonicalMap',canonicalMap, ...
    'method',method,'classification','circle_spectrum_only');
end
