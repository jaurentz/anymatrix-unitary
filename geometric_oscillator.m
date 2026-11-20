function [U,properties] = geometric_oscillator(n,h,omega0,kappa,method,mode)
%GEOMETRIC_OSCILLATOR Orthogonal-symplectic steps for a periodic oscillator chain.
%   U=GEOMETRIC_OSCILLATOR(N,H,OMEGA0,KAPPA,METHOD) returns the real
%   2N-by-2N orthogonal-symplectic representative of an integration step
%   for q''+K*q=0, with K=OMEGA0^2*I+KAPPA*(2*I-P-P').
%   P is the cyclic shift, OMEGA0>0, KAPPA>=0, H>0. Defaults:
%   H=0.5, OMEGA0=1, KAPPA=1, METHOD='midpoint'.
%   METHOD='midpoint', 'gauss2', 'gauss3' (s-stage Gauss-Legendre)
%   or 'verlet' (velocity/Stoermer-Verlet, stable H*max(omega)<2).
%   METHOD='euler_kick_drift' or 'euler_drift_kick' are the
%   two distinct first-order symplectic-Euler physical updates.
%   They preserve a positive quadratic energy with cross-block terms
%   when H*max(omega)<2; their normalized orthogonal representative
%   has the same exact phase multiset as Verlet, but not the same M,G.
%
%   The real physical step M is symplectic and preserves a positive
%   quadratic form G: M'*G*M=G. The returned U=T*M/T is also
%   orthogonal and symplectic, where G=T'*T. T is not required to be
%   symplectic. The eigenvalues of U and of the physical step are the
%   exact conjugate pairs [z;conj(z)], for frequencies omega_j:
%       omega_j^2=OMEGA0^2+4*KAPPA*sin(pi*j/N)^2.
%   The two Euler orderings and Verlet share the phase pairs
%   exp(+/-2i*asin(H*omega_j/2)), but have different phase-space steps.
%   Gauss: z=R_s(-1i*H*omega_j), R_s=P_s(z)/P_s(-z), the diagonal
%   Pade approximant to exp(z). Verlet:
%       z=1-(H*omega_j)^2/2 -1i*H*omega_j*sqrt(1-(H*omega_j)^2/4).
%
%   Terminal 'operator' returns O(X,FLAG) with O(n log n) FFT actions.
%   Terminal 'pencil' returns struct with A=G*M, B=G, M, G, T, U.
%   A-lambda*B has unit-circle eigenvalues but is NOT generally a
%   strict Euclidean-unitary pencil: B\A=M need not be orthogonal.
%   Terminal 'polynomial' (Verlet only) returns symmetric palindromic
%   P(lambda)=lambda^2*I+lambda*(H^2*K-2I)+I (coefficients A0,A1,A2).
%   Terminal 'spectral' returns the standard finite trace-spectrum descriptor.
%
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration,
%   2nd ed., 2006. [97] https://doi.org/10.1007/3-540-30666-8
%   E. Hairer, C. Lubich, G. Wanner, Acta Numerica 12 (2003), 399-450.
%   [98] https://doi.org/10.1017/S0962492902000144
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0, U=[]; return; end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(h), h=0.5; end
if nargin<3 || isempty(omega0), omega0=1; end
if nargin<4 || isempty(kappa), kappa=1; end
if nargin<5 || isempty(method), method='midpoint'; end
if nargin<6 || isempty(mode), mode='matrix'; end
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(omega0,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(kappa,{'numeric'},{'scalar','real','finite','nonnegative'});
if ~(ischar(method)||isstring(method)) || ...
        ~(ischar(mode)||isstring(mode))
    error('geometric_oscillator:Mode','METHOD and MODE must be text.');
end
method=lower(char(method)); mode=lower(char(mode));
if ~any(strcmp(method,{'midpoint','gauss2','gauss3','verlet', ...
        'euler_kick_drift','euler_drift_kick'}))
    error('geometric_oscillator:Method','Unknown integration method.');
end
if ~any(strcmp(mode,{'matrix','operator','pencil','polynomial','spectral'}))
    error('geometric_oscillator:Mode','Unknown output mode.');
end
if strcmp(mode,'polynomial') && ~strcmp(method,'verlet')
    error('geometric_oscillator:Polynomial', ...
        'The implemented palindromic position recurrence is for Verlet.');
end
j=(0:n-1)';
w=sqrt(omega0^2+4*kappa*sin(pi*j/n).^2);
if any(strcmp(method,{'verlet','euler_kick_drift','euler_drift_kick'}))
    if h*max(w)>=2
        error('geometric_oscillator:Stability', ...
            'Verlet and symplectic Euler require h*max(omega_j)<2 for a positive metric.');
    end
    x=h*w;
    z=(1-0.5*x.^2)-1i*x.*sqrt(1-0.25*x.^2);
else
    switch method
        case 'midpoint', coeff=[1 1/2];
        case 'gauss2',   coeff=[1 1/2 1/12];
        case 'gauss3',   coeff=[1 1/2 1/10 1/120];
    end
    z=polyval(fliplr(coeff),-1i*h*w) ./ ...
      polyval(fliplr(coeff), 1i*h*w);
end
if strcmp(mode,'operator')
    U=@(X,flag)um_geometric_oscillator_apply(X,flag,z);
    return
end
if strcmp(mode,'spectral')
    U=um_geometric_spectral('geometric_oscillator',z,w);
    return
end
P=sparse(1:n,[2:n,1],1,n,n);
K=omega0^2*speye(n)+kappa*(2*speye(n)-P-P');
if strcmp(mode,'polynomial')
    I=eye(n);
    U=struct('A0',I,'A1',h*h*full(K)-2*I,'A2',I, ...
        'stiffness',K,'dimension',n,'method',method);
    return
end
% The frequency sequence is conjugation-symmetric in Fourier space.
% Its inverse FFT is the first column of a real symmetric circulant.
ind=mod((0:n-1)'-(0:n-1),n)+1;
first=ifft(z);
C=real(first(ind)); S=-imag(first(ind));
U=[C,S;-S,C];
if strcmp(mode,'matrix'), return; end
% In canonical phase-space coordinates x=[q;p], G=diag(Gq,I).
if any(strcmp(method,{'euler_kick_drift','euler_drift_kick'}))
    % For symplectic Euler, G has genuine q-p cross terms. The
    % symmetric positive square root gives the canonical oriented
    % 2-by-2 phase rotation in each Fourier stiffness mode.
    I=eye(n); K=full(K);
    if strcmp(method,'euler_kick_drift')
        % Kick p_+=p-hKq, then drift q_+=q+h*p_+.
        M=[I-h*h*K,h*I;-h*K,I];
        C=-0.5*h*K;
    else
        % Drift q_+=q+hp, then kick p_+=p-hKq_+.
        M=[I,h*I;-h*K,I-h*h*K];
        C=0.5*h*K;
    end
    G=[K,C;C,I];
    [V,D]=eig(G);
    diagonal=real(diag(D));
    if any(diagonal<=0)
        error('geometric_oscillator:Metric','Expected a positive Euler metric.');
    end
    T=V*diag(sqrt(diagonal))*V';
    T=real((T+T')/2);  % symmetric positive square root
else
    if strcmp(method,'verlet')
        Gq=K*(speye(n)-0.25*h*h*K);  % modified energy
        D=full(speye(n)-0.5*h*h*K);
        M=[D,h*eye(n);-h*full(Gq),D];
        tcoeff=w.*sqrt(1-0.25*(h*w).^2);
    else
        Gq=K;                         % original physical energy
        tcoeff=w;
    end
    Tq=real(ifft(tcoeff));
    T=blkdiag(Tq(ind),eye(n));        % G=T'*T (up to rounding)
    G=blkdiag(full(Gq),eye(n));
    if ~strcmp(method,'verlet')
        M=T\(U*T);                    % the Gauss physical step in (q,p)
    end
end
A=G*M; B=G;
U=struct('A',A,'B',B,'M',M,'G',G,'T',T,'U',U, ...
    'classification','circle_spectrum_only','method',method);
end
