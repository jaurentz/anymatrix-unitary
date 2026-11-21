function [U,properties]=hermite_oscillator(n,theta,offset,mode)
%HERMITE_OSCILLATOR Finite Gauss-Hermite harmonic-oscillator unitary flow.
%   U=HERMITE_OSCILLATOR(N,THETA,OFFSET,MODE) returns the N-by-N
%   strictly unitary sample-space model U(theta)=Q*diag(
%   exp(-1i*theta*(k+OFFSET)))*Q', k=0,...,N-1, where
%   Q=GAUSS_HERMITE_TRANSFORM(N) is the real orthogonal Gauss-Hermite
%   QUADRATURE EVALUATION matrix with rows indexed by increasing
%   Gauss-Hermite nodes x_j and columns by normalized physicists'
%   Hermite-polynomial degree k. Q(j,k+1)=sqrt(w_j)*p_k(x_j),
%   where p_k are orthonormal for exp(-x^2)/sqrt(pi) and w_j are
%   quadrature weights normalized to sum 1.
%
%   This is the strict finite oscillator spectral model in quadrature
%   sample coordinates, NOT the Candan-Kutay-Ozaktas DFT-compatible
%   fractional Fourier transform, NOT the asymptotic FFT-based algorithm
%   of Campos et al., and NOT a claim U(pi/2)=DFT(N).
%   The normalized trace eigenphases are EXACTLY
%   { exp(-1i*theta*(k+OFFSET)):k=0,...,N-1 }.
%
%   OFFSET=0 (default) omits the global ground-state phase, giving
%   U(theta+2*pi)=U(theta) and U(0)=I. OFFSET=1/2 includes
%   the quantum harmonic-oscillator zero-point energy, producing
%   U(theta+2*pi)=-U(theta) and U(theta+4*pi)=U(theta).
%   Both choices obey U(theta)*U(phi)=U(theta+phi) and adjoint inversion.
%   N=1..256, THETA real scalar (default .7), OFFSET=0 or 1/2.
%
%   MODE='matrix'      default dense n-by-n unitary.
%   MODE='spectral'    exact finite phase multiset and moments
%                       without numerical eig(U), no asserted limit.
%   MODE='hamiltonian' real symmetric n-by-n discrete number
%                       operator H=Q*diag(k+OFFSET)*Q';
%                       expm(-1i*theta*H)=U(theta).
%   MODE='basis'       real Gauss-Hermite evaluation Q, quadrature
%                       nodes, normalized positive weights, H spectrum.
%
%   Construction cost is dense O(N^3) eigenvector setup and O(N^3)
%   multiplication; no unwarranted matrix-free/FFT speedup is claimed.
%   Differentiation: dU/dtheta=-1i*H*U.
%
%   [104] Golub-Welsch, Math. Comp. 23 (1969), quadrature eigenbasis.
%         https://doi.org/10.1090/S0025-5718-69-99647-1
%   [143] R. G. Campos, L. Z. Juarez, J. S. Dehesa (2012), SIAM JSC,
%         a DISTINCT fast fractional Fourier quadrature algorithm.
%         https://doi.org/10.1137/100812677
%   [144] Candan, Kutay and Ozaktas (2000), IEEE TSP, a DISTINCT
%         DFT-compatible discrete fractional Fourier definition.
%         https://doi.org/10.1109/78.839980
properties={'square','unitary','normal','complex', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
        &&n==floor(n)&&n>=1&&n<=256)
    error('hermite_oscillator:Order','N must be 1..256.');
end
if nargin<2||isempty(theta),theta=.7;end
if nargin<3||isempty(offset),offset=0;end
if nargin<4||isempty(mode),mode='matrix';end
if ~(isnumeric(theta)&&isscalar(theta)&&isreal(theta) ...
       &&isfinite(theta))
    error('hermite_oscillator:Theta','THETA must be real finite scalar.');
end
if ~(isnumeric(offset)&&isscalar(offset)&&isreal(offset) ...
       &&any(offset==[0,.5]))
    error('hermite_oscillator:Offset','OFFSET must be 0 or 1/2.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('hermite_oscillator:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','spectral','basis','hamiltonian'}))
    error('hermite_oscillator:Mode','Unknown output mode.');
end
energy=(0:n-1)'+offset;
lambda=exp(-1i*theta*energy);
if strcmp(mode,'spectral')
    blank=struct('kind','unknown','density',[],'support',[], ...
       'atoms',[],'moments',[]);
    S=struct('family','hermite_oscillator','dimension',n, ...
       'normalization','probability (density with respect to dtheta)', ...
       'finite',struct('kind','exact','eigenvalues',@()lambda, ...
        'moments',@(k)hermite_moments(lambda,k), ...
        'atoms',[]), ...
       'limit',blank,'ensemble',blank,'underlying',blank, ...
       'reference','Golub-Welsch quadrature and finite spectral harmonic-oscillator model', ...
       'notes',['Exact finite oscillator energy eigenvalues k+OFFSET ', ...
          'and phase eigenvalues exp(-i theta(k+OFFSET)). ', ...
          'No equality at theta=pi/2 with normalized DFT is claimed: ', ...
          'Gauss-Hermite quadrature eigenvectors differ from the ', ...
          'DFT-compatible discrete Hermite-like eigenvectors of ', ...
          'Candan, Kutay and Ozaktas. No asymptotic or ensemble ', ...
          'trace density is asserted.']);
    S.oscillator_energies=energy;
    S.theta=theta;S.offset=offset;
    U=S;
    return
end
Q=gauss_hermite_transform(n);
switch mode
    case 'matrix'
        U=(Q.*lambda.')*Q';
    case 'hamiltonian'
        U=(Q.*energy.')*Q';
        U=real((U+U')/2);
    case 'basis'
        diagonal=zeros(n,1);
        offdiag=sqrt((1:n-1)'/2);
        Jacobi=diag(diagonal)+diag(offdiag,1)+diag(offdiag,-1);
        nodes=diag(Q*Jacobi*Q');
        U=struct('Q',Q,'nodes',nodes, ...
          'weights',Q(:,1).^2, ...
          'oscillator_energies',energy, ...
          'quadrature_measure','exp(-x^2) dx / sqrt(pi)', ...
          'basis_polynomials','orthonormal physicists Hermite', ...
          'representation','quadrature-weighted samples by increasing node');
end
end

function moments=hermite_moments(lambda,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
shape=size(k);
values=reshape(k,1,[]);
moments=reshape(mean(bsxfun(@power,lambda(:),values),1),shape);
end
