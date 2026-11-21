function [U,properties]=discrete_fractional_fourier(n,alpha,mode)
%DISCRETE_FRACTIONAL_FOURIER Candan--Kutay--Ozaktas exact discrete FrFT.
%   U=DISCRETE_FRACTIONAL_FOURIER(N,ALPHA) returns the N-by-N unitary
%   discrete fractional Fourier transform of Candan, Kutay and Ozaktas
%   (IEEE TSP 48 (2000), Section III, equations (9), (15), (22)).
%   Its N-by-N normalized DFT uses F(j,k)=exp(-2*pi*i*j*k/N)/sqrt(N).
%   An orthonormal REAL discrete Hermite-like basis V simultaneously
%   diagonalizes the DFT F and the real symmetric discrete oscillator
%       S = periodic_shift + periodic_shift' +
%           diag(2*cos(2*pi*j/N)-4),  j=0,...,N-1.
%   S has a doubly degenerate zero eigenvalue when N is divisible by 4.
%   The basis is therefore obtained from *separate even and odd parity
%   blocks* (not from a generic eig(S) call). Descending eigenvalue
%   ordering in each parity block gives the zero-crossing index:
%     N odd:  k = 0,1,...,N-1;
%     N even: k = 0,1,...,N-2,N (N-1 is SKIPPED).
%   This exceptional even-N index is essential for F_1=DFT(N).
%
%   The transform is U = V*diag(exp(-i*pi*alpha*k/2))*V'.
%   Exact-mathematics identities: U(0)=I, U(1)=DFT(N),
%   U(alpha)*U(beta)=U(alpha+beta), U(alpha)'=U(-alpha),
%   U(alpha+4)=U(alpha), and strict unitarity.
%   Basis signs choose largest-magnitude sample positive; this does
%   not change U. This is *not* the Gauss-Hermite oscillator transform.
%
%   MODE='matrix' (default): N-by-N explicit unitary.
%        'basis': struct V, oscillator S, discrete Hermite indices,
%                 parity, signed commuting-block eigenvalues.
%        'oscillator': real symmetric discrete difference S (NOT unitary).
%        'hamiltonian': real symmetric number operator V*diag(k)*V',
%                       such that U=expm(-i*pi*alpha*H/2).
%        'spectral': EXACT finite phase multiset and trace moments.
%   N<=256 to bound dense O(N^3) eigensolver construction. No 'operator'
%   interface is offered because a cached dense basis would *not* be
%   genuinely faster or subquadratic storage. This is the spectral
%   Candan discrete FrFT, not an O(N log N) approximation to the
%   continuous fractional transform integral.
%
%   [161] C. Candan, M. A. Kutay, H. M. Ozaktas,
%   "The discrete fractional Fourier transform", IEEE TSP 48 (2000).
%   https://doi.org/10.1109/78.839980
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',256});
if nargin<2||isempty(alpha),alpha=.37;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
if nargin<3||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('discrete_fractional_fourier:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','basis','oscillator','hamiltonian','spectral'}))
    error('discrete_fractional_fourier:Mode', ...
        'Use matrix, basis, oscillator, hamiltonian, or spectral.');
end
if strcmp(mode,'spectral')
    k=hermite_indices(n);
    phases=exp(-.5i*pi*alpha*k);
    U=unitary_spectral_reference('discrete_fractional_fourier',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()phases;
    U.finite.moments=@(orders)phase_moments(phases,orders);
    U.finite.atoms=struct('angles',angle(phases), ...
        'weights',ones(n,1)/n);
    U.hermite_indices=k;
    U.order=alpha;
    U.reference='Candan--Kutay--Ozaktas, IEEE TSP 48 (2000), Sec. III';
    U.notes=['The DFT-compatible discrete Hermite branch skips index ', ...
        'N-1 for even N and uses index N instead. Exact finite trace ', ...
        'phases are exp(-i*pi*alpha*k/2). This spectral expansion ', ...
        'is distinct from continuous-integral quadrature approximations.'];
    return
end
S=discrete_oscillator(n);
if strcmp(mode,'oscillator'),U=S;return;end
[V,k,oscillatorEigenvalues,parity]=joint_basis(S);
if strcmp(mode,'basis')
    U=struct('V',V,'hermite_indices',k,'parity',parity, ...
        'commuting_oscillator',S, ...
        'oscillator_eigenvalues',oscillatorEigenvalues, ...
        'normalization','real orthonormal vectors in natural DFT sample ordering', ...
        'ordering','descending oscillator eigenvalues separately by parity');
    return
end
if strcmp(mode,'hamiltonian')
    U=real((V.*k.')*V');
    U=(U+U')/2;
    return
end
lambda=exp(-.5i*pi*alpha*k);
U=(V.*lambda.')*V';
end

function S=discrete_oscillator(n)
j=(0:n-1)';
S=diag(2*cos(2*pi*j/n)-4);
for m=1:n
    S(m,mod(m,n)+1)=S(m,mod(m,n)+1)+1;
    S(m,mod(m-2,n)+1)=S(m,mod(m-2,n)+1)+1;
end
S=(S+S')/2;
end

function k=hermite_indices(n)
if mod(n,2)
    k=(0:n-1)';
else
    k=[(0:n-2)';n];
end
end

function [V,k,eigvalues,parity]=joint_basis(S)
n=size(S,1);
% Construct orthonormal even and odd parity subspaces of f(j)=f(-j).
% For even n the sample at n/2 is a second fixed point.
Ve=zeros(n,floor(n/2)+1);
Ve(1,1)=1;
q=1;
for j=1:floor((n-1)/2)
    q=q+1;
    Ve(j+1,q)=1/sqrt(2);
    Ve(n-j+1,q)=1/sqrt(2);
end
if mod(n,2)==0
    Ve(n/2+1,end)=1;
end
Vo=zeros(n,floor((n-1)/2));
for j=1:size(Vo,2)
    Vo(j+1,j)=1/sqrt(2);
    Vo(n-j+1,j)=-1/sqrt(2);
end
% In the odd-n case Ve has the right dimension floor(n/2)+1;
% for even n it has n/2+1, including the Nyquist parity fixed point.
V=zeros(n);eigvalues=zeros(n,1);parity=zeros(n,1);
jcol=0;
for which=1:2
    if which==1, P=Ve;signParity=1;else,P=Vo;signParity=-1;end
    count=size(P,2);
    if count==0,continue;end
    T=P'*S*P;
    T=(T+T')/2;
    [Z,D]=eig(T);
    [lambda,perm]=sort(real(diag(D)),'descend');
    Z=real(Z(:,perm));
    W=P*Z;
    for j=1:count
        [~,pivot]=max(abs(W(:,j)));
        if W(pivot,j)<0,W(:,j)=-W(:,j);end
        if signParity==1,idx=2*(j-1);else,idx=2*j-1;end
        target=find(hermite_indices(n)==idx,1);
        if isempty(target)
            error('discrete_fractional_fourier:Index','Invalid Hermite index.');
        end
        V(:,target)=W(:,j);
        eigvalues(target)=lambda(j);
        parity(target)=signParity;
        jcol=jcol+1;
    end
end
if jcol~=n
    error('discrete_fractional_fourier:Basis','Incomplete parity basis.');
end
k=hermite_indices(n);
end

function moments=phase_moments(z,orders)
validateattributes(orders,{'numeric'},{'real','integer','finite'});
shape=size(orders);p=reshape(orders,1,[]);
moments=reshape(mean(bsxfun(@power,z(:),p),1),shape);
end
