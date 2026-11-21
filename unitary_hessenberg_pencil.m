function [U,properties]=unitary_hessenberg_pencil(n,poles,seed,variant,mode)
%UNITARY_HESSENBERG_PENCIL Core-factored strict unitary Hessenberg pair.
%   U=UNITARY_HESSENBERG_PENCIL(N,POLES,SEED,VARIANT,MODE) constructs
%   independent unitary upper-Hessenberg A,B and returns U=B'*A.
%   A and B are individually unitary, so A-z*B is a strict
%   unitary-unitary generalized eigenvalue pencil:
%       eig(A,B)=eig(B'*A),  |eig(B'*A)|=1.
%   This is not the larger (nonunitary) Hessenberg--unitary-
%   Hessenberg pencil from generic rational QZ/pole swapping.
%
%   POLES has N-1 projective ratios
%       xi_j = A(j+1,j)/B(j+1,j), j=1,...,N-1.
%   For each prescribed finite nonzero pole with modulus s,
%   choose core sines rhoA=t*min(1,s), rhoB=t*min(1,1/s),
%   t=0.7, with strict-disk core moduli sqrt(1-rho^2).
%   The left row phases of A encode arg(xi_j). Independent Schur
%   phases in both core products and terminal phases are seeded.
%   A=diag(leftA)*H(alphaA,betaA), B=H(alphaB,betaB), where H
%   is the unitary upper-Hessenberg ordered product of 2-by-2
%   Verblunsky cores. Compact parameters require O(N) storage.
%
%   A finite zero or Inf pole sets respectively rhoA=0 or
%   rhoB=0: the appropriate core is exactly diagonal, so the
%   associated Hessenberg factor is REDUCIBLE. They are legitimate
%   projective ratios (0/nonzero and nonzero/0), but not proper
%   unreduced Hessenberg pairs. They cannot be obtained with all
%   interior |alpha_j|<1. Finite nonzero poles are limited to
%   1e-4<=abs(xi)<=1e4 to resolve subdiagonal ratios in double
%   precision; a final parameters diagnostic provides achieved
%   ratios including roundoff.
%
%   VARIANT='complex' (default) permits arbitrary complex poles;
%   VARIANT='real' requires real finite poles or Inf and returns
%   real orthogonal A,B,U. Independent real signs are seeded.
%   SEED (default 1) preserves the caller's RNG state.
%   Default POLES=linspace(.4,2.5,N-1)' is proper and real.
%
%   MODE='matrix' (default): explicit unitary quotient U.
%        'pencil': struct A, B, U, projective poles, proper flags.
%        'parameters': O(N) core factors, closure phases, pole data.
%        'operator': O(N) forward and adjoint U(X,flag), multi-RHS,
%                    without constructing A, B or U.
%        'spectral': finite trace eigenphase law unknown in general;
%                    poles are NOT eigenvalues of the pencil.
%   Explicit matrices are capped at N<=512; 'operator' and
%   'parameters' support N<=2e6 and use O(N) factor data.
%
%   [163] Camps, Mach, Vandebril, Watkins, Pole-Swapping Algorithms
%   for the Eigenvalue Problem, SIAM Spotlights (2025), Chapters 4-5.
%   https://doi.org/10.1137/1.9781611978377
%   [164] Camps, Mach, Vandebril, Watkins, ETNA 52 (2020), 480-508.
%   https://doi.org/10.1553/etna_vol52s480
%   The explicit strict pair is a derived core-product specialization,
%   not a claim that either reference implemented this generator.
properties={'square','unitary','normal','scalable', ...
    'parameter-dependent','random'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'}, ...
    {'scalar','integer','finite','>=',2,'<=',2e6});
if nargin<2||isempty(poles),poles=linspace(.4,2.5,n-1).';end
if nargin<3||isempty(seed),seed=1;end
if nargin<4||isempty(variant),variant='complex';end
if nargin<5||isempty(mode),mode='matrix';end
validateattributes(seed,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
if ~(ischar(variant)||(isstring(variant)&&isscalar(variant))) ...
        || ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('unitary_hessenberg_pencil:Mode', ...
        'VARIANT and MODE must be scalar text.');
end
variant=lower(char(variant));mode=lower(char(mode));
if ~any(strcmp(variant,{'real','complex'}))
    error('unitary_hessenberg_pencil:Variant','Unknown variant.');
end
if ~any(strcmp(mode,{'matrix','pencil','parameters','operator','spectral'}))
    error('unitary_hessenberg_pencil:Mode','Unknown output mode.');
end
if ~isnumeric(poles)||~isvector(poles)||numel(poles)~=n-1 ...
        ||any(isnan(poles(:)))|| ...
        any(~(isfinite(poles(:))|isinf(poles(:))))|| ...
        any(isinf(poles(:))&imag(poles(:))~=0)
    error('unitary_hessenberg_pencil:Poles', ...
        'N-1 finite poles or real Inf must be supplied.');
end
poles=poles(:);
if strcmp(variant,'real') && any(abs(imag(poles(isfinite(poles))))>1e-13)
    error('unitary_hessenberg_pencil:Real', ...
        'Real orthogonal factors require real projective poles.');
end
mag=abs(poles);mask=mag>0&isfinite(mag);
if any(mag(mask)<1e-4|mag(mask)>1e4)
    error('unitary_hessenberg_pencil:ResolvedPoles', ...
        'Finite nonzero pole moduli must lie in [1e-4,1e4].');
end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('unitary_hessenberg_pencil',n);
    U.poles=poles;
    U.reference='Camps--Mach--Vandebril--Watkins, ETNA 52 (2020)';
    U.notes=['Strict unitary--unitary upper-Hessenberg pencil A-zB. ', ...
        'Pole ratios are subdiagonal projective quantities, not ', ...
        'eigenvalues/eigenphases. Generic finite trace law unknown.'];
    return
end
if any(strcmp(mode,{'matrix','pencil'}))&&n>512
    error('unitary_hessenberg_pencil:ExplicitSize', ...
        'Dense output is capped at order 512; use operator/parameters.');
end
prior=rng;restore=onCleanup(@()rng(prior)); %#ok<NASGU>
rng(seed,'twister');
t=.7;
rhoA=t*min(ones(n-1,1),mag);
rhoB=t*min(ones(n-1,1),1./mag);
rhoB(isinf(mag))=0;
ampA=sqrt((1-rhoA).*(1+rhoA));
ampB=sqrt((1-rhoB).*(1+rhoB));
if strcmp(variant,'real')
    sgna=2*(rand(n-1,1)>.5)-1;
    sgnb=2*(rand(n-1,1)>.5)-1;
    alphaA=ampA.*sgna;
    alphaB=ampB.*sgnb;
    betaA=2*(rand>.5)-1;
    betaB=2*(rand>.5)-1;
    phaseA=ones(n,1);
    p=real(poles);
    phaseA(2:end)=sign(p);
    phaseA(2:end)=phaseA(2:end)+(phaseA(2:end)==0);
else
    alphaA=ampA.*exp(2i*pi*rand(n-1,1));
    alphaB=ampB.*exp(2i*pi*rand(n-1,1));
    betaA=exp(2i*pi*rand);
    betaB=exp(2i*pi*rand);
    phaseA=[1;exp(1i*angle(poles))];
end
% Diagonal reducible cores need exact unimodular alpha; using an
% arbitrary rounded exp(i*theta) could produce spurious O(sqrt(eps))
% subdiagonals and destroy exact zero/infinite projective poles.
alphaA(rhoA==0)=1;
alphaB(rhoB==0)=1;
if any(abs(alphaA(rhoA>0))>=1)||any(abs(alphaB(rhoB>0))>=1)
    error('unitary_hessenberg_pencil:Roundoff', ...
        'Core parameter too close to unit modulus to resolve.');
end
rhoAeff=sqrt((1-abs(alphaA)).*(1+abs(alphaA)));
rhoBeff=sqrt((1-abs(alphaB)).*(1+abs(alphaB)));
achieved=phaseA(2:end).*(rhoAeff./rhoBeff);
achieved(rhoBeff==0)=Inf;
achieved(rhoAeff==0 & rhoBeff>0)=0;
proper=(rhoAeff>0 & rhoBeff>0);
if strcmp(mode,'parameters')
    U=struct('alphaA',alphaA,'alphaB',alphaB, ...
        'betaA',betaA,'betaB',betaB,'left_phases_A',phaseA, ...
        'target_poles',poles,'achieved_poles',achieved, ...
        'rhoA',rhoAeff,'rhoB',rhoBeff,'proper',proper, ...
        'is_proper',all(proper),'dimension',n,'seed',seed, ...
        'variant',variant,'core_count',2*(n-1));
    return
end
[~,steps]=um_core_order(n,'hessenberg',[]);
if strcmp(mode,'operator')
    U=@(X,flag)um_hessenberg_pencil_apply(X,flag, ...
        alphaA,betaA,phaseA,alphaB,betaB,steps);
    return
end
A=bsxfun(@times,phaseA, ...
    um_core_apply(eye(n),'notransp',alphaA,betaA,steps));
B=um_core_apply(eye(n),'notransp',alphaB,betaB,steps);
T=B'*A;
if strcmp(mode,'matrix'),U=T;return;end
U=struct('A',A,'B',B,'U',T, ...
    'poles',achieved,'target_poles',poles,'proper',proper, ...
    'is_proper',all(proper),'dimension',n, ...
    'note',['A and B individually unitary upper Hessenberg; ', ...
        'projective zero/Inf subdiagonal poles imply reducibility.']);
end
