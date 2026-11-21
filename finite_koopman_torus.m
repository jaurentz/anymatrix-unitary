function [U,properties]=finite_koopman_torus(m,M,mode,translation)
%FINITE_KOOPMAN_TORUS Exact finite Koopman pullback on a 2D modular torus.
%   U=FINITE_KOOPMAN_TORUS(M,MATRIX,MODE,TRANSLATION) is a real
%   orthogonal N-by-N permutation (N=M^2), representing the pullback
%       (U*f)(x)=f(T(x)),   T(x)=MATRIX*x+TRANSLATION mod M,
%   on (Z/M Z)^2 with normalized counting measure. The default MATRIX
%       [2 1; 1 1]
%   is the finite classical Arnold cat automorphism (det=1), NOT the
%   metaplectic QUANTUM_CAT_MAP. TRANSLATION defaults to [0;0].
%   Alternative MATRIX [1 1;0 1] yields a shear; MATRIX=eye(2)
%   with nonzero TRANSLATION yields finite torus translations.
%
%   Grid label x=(x1,x2), x1,x2=0,...,M-1, is stored as
%       i=1+x1+M*x2.
%   Input MATRIX must be integer 2-by-2 with determinant coprime to M,
%   ensuring bijectivity of the affine T. This is checked explicitly.
%   The finite Koopman is unitary because its states are permuted,
%   not because an arbitrary finite projection of an infinite
%   Koopman operator remains unitary.
%
%   MODE='matrix' (default): sparse permutation, N<=4096.
%   MODE='operator': O(N) permutation forward/adjoint, N<=1e6.
%   MODE='map': map indices, inverse indices, cycles, lengths, fixed
%       points and dual-frequency map for TRANSLATION=0.
%   MODE='spectral': EXACT finite cycle roots and normalized trace
%       moments tr(U^k)/N = # fixed points of T^k / N for integer k.
%       Generic continuous-torus Koopman spectral types and their
%       limiting laws are NOT implied by this finite calculation.
%
%   With TRANSLATION=0, Fourier characters
%       chi_k(x)=exp(2*pi*1i*k'*x/M)
%   obey U chi_k=chi_(MATRIX'*k). The normalized full 2D DFT
%   transforms U into a *dual permutation* with that convention.
%   For a translated map, an additional phase exp(2*pi*i*k'*b/M)
%   appears in the dual action.
%
%   [93] Hannay and Berry, Quantization of Linear Maps on a Torus,
%       Physica D 1 (1980), 267-290.
%       https://doi.org/10.1016/0167-2789(80)90026-3
%   [146] Colbrook and Townsend, rigorous Koopman spectra, CPAM 2024.
%       https://doi.org/10.1002/cpa.22125
%   This particular exact finite permutation family is elementary
%   modular dynamics; it is not claimed as a published example of [146].
properties={'square','unitary','normal','real','orthogonal', ...
            'sparse','scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(m)&&isscalar(m)&&isreal(m)&&isfinite(m) ...
        &&m==fix(m)&&m>=2&&m<=1000)
    error('finite_koopman_torus:Order','M must be an integer 2..1000.');
end
if nargin<2||isempty(M),M=[2 1;1 1];end
if ~(isnumeric(M)&&isreal(M)&&isequal(size(M),[2 2]) ...
        &&all(isfinite(M(:)))&&all(M(:)==fix(M(:))) ...
        &&max(abs(M(:)))<=1e6)
    error('finite_koopman_torus:Matrix','Matrix must be a finite integer 2-by-2 with entries <= 1e6.');
end
if nargin<3||isempty(mode),mode='matrix';end
if nargin<4||isempty(translation),translation=[0;0];end
if ~(isnumeric(translation)&&isreal(translation) ...
        &&isvector(translation)&&numel(translation)==2 ...
        &&all(isfinite(translation(:))) ...
        &&all(translation(:)==fix(translation(:))) ...
        &&max(abs(translation(:)))<=1e6)
    error('finite_koopman_torus:Translation','Translation must be an integer two-vector.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('finite_koopman_torus:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','map','spectral'}))
    error('finite_koopman_torus:Mode','Unsupported output mode.');
end
detM=double(M(1,1))*double(M(2,2)) ...
     -double(M(1,2))*double(M(2,1));
if gcd(abs(detM),double(m))~=1
    error('finite_koopman_torus:Invertibility', ...
          'Determinant must be coprime to M for a bijective map.');
end
M=mod(double(M),double(m));
b=mod(double(translation(:)),double(m));
N=double(m)*double(m);
if strcmp(mode,'matrix')&&N>4096
    error('finite_koopman_torus:ExplicitOrder', ...
        'Explicit sparse matrix supports M^2<=4096; use operator mode.');
end
% Exactly representable integer products: M reduced modulo m, m<=1000.
j=(0:N-1)';
x=mod(j,m);y=floor(j/m);
X=mod(M(1,1)*x+M(1,2)*y+b(1),m);
Y=mod(M(2,1)*x+M(2,2)*y+b(2),m);
perm=1+X+m*Y;
% A single bijective map index is used throughout matrix, operator,
% and spectral modes; never cast a noninvertible map as unitary.
invperm=zeros(N,1);invperm(perm)=1:N;
if any(invperm==0)
    error('finite_koopman_torus:Invertibility','Map has repeated states.');
end
if strcmp(mode,'matrix')
    U=sparse((1:N)',perm,1,N,N);
    return
end
if strcmp(mode,'operator')
    U=@(A,flag)um_finite_koopman_apply(A,flag,perm,invperm);
    return
end
[cycles,lengths]=um_finite_koopman_cycles(perm);
if strcmp(mode,'map')
    % k -> M^T*k dual-frequency permutation in character coordinates.
    K1=mod(M(1,1)*x+M(2,1)*y,m);
    K2=mod(M(1,2)*x+M(2,2)*y,m);
    dual=1+K1+m*K2;
    dualPhase=exp(2i*pi*(x*b(1)+y*b(2))/m);
    U=struct('dimension',N,'modulus',m,'matrix',M, ...
      'translation',b,'permutation',perm,'inverse',invperm, ...
      'cycles',{cycles},'cycle_lengths',lengths, ...
      'dual_frequency_permutation',dual, ...
      'dual_frequency_phases',dualPhase, ...
      'fixed_points',sum(perm==(1:N)'), ...
      'order',double(lcm_of_lengths(lengths)), ...
      'convention','(U f)(x)=f(T(x)), x1 fastest; characters exp(+2pi i k*x/m)');
    return
end
blank=struct('kind','unknown','density',[], ...
        'support',[],'atoms',[],'moments',[]);
U=struct('family','finite_koopman_torus','dimension',N, ...
      'normalization','probability (density with respect to dtheta)', ...
      'finite',struct('kind','exact', ...
        'eigenvalues',@()um_finite_koopman_cycle_eigs(lengths), ...
        'moments',@(k)um_finite_koopman_cycle_moments(lengths,N,k), ...
        'atoms',[]), ...
      'limit',blank,'ensemble',blank,'underlying',blank, ...
      'reference','Finite modular Koopman pullback; Colbrook and Townsend (2024) for infinite operator distinction', ...
      'cycle_lengths',lengths,'fixed_points',sum(perm==(1:N)'), ...
      'map_order',double(lcm_of_lengths(lengths)), ...
      'notes',['Exact finite permutation roots and normalized trace ', ...
       'moments, with all cycle multiplicities. Spectrum belongs to ', ...
       'roots of unity and differs from the infinite Koopman operator ', ...
       'on a continuous chaotic torus. This is NOT metaplectic ', ...
       'quantization of a cat map.']);
end

function v=lcm_of_lengths(lengths)
% Avoid silent integer precision loss for possibly very large order.
v=uint64(1);
for j=1:numel(lengths)
    t=uint64(lengths(j));
    divisor=gcd(v,t);
    if double(v)/double(divisor)>double(intmax('uint64'))/double(t)
        v=uint64(0);return
    end
    v=(v/divisor)*t;
end
end
