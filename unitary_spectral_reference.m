function S = unitary_spectral_reference(family,varargin)
%UNITARY_SPECTRAL_REFERENCE Analytic spectral reference data for unitary generators.
%   S=UNITARY_SPECTRAL_REFERENCE(FAMILY,PARAMETERS...) returns a standardized
%   descriptor, without constructing the full matrix or advancing the caller's RNG state.
%   Seeded Haar tensor families generate only small base-factor matrices
%   and compute their eigenvalues, then derive large-matrix spectra analytically.
%   Public generators expose this through their terminal 'spectral' mode:
%      S=anymatrix('unitary/dirac_automaton',N,M,'spectral');
%
%   S.FINITE describes the equally-weighted finite trace measure.
%   S.LIMIT describes the limiting equally-weighted trace measure, if known.
%   S.ENSEMBLE describes an ensemble expectation (not a sample).
%   S.UNDERLYING describes an infinite OPUC distinguished-vector measure,
%   which is NOT generally a finite-matrix trace measure.
%   Densities are relative to dtheta, not dtheta/(2*pi).
%   All phase supports use [-pi,pi]. Atoms have fields angles and weights.
%   FINITE.EIGENVALUES is a zero-argument handle; MOMENTS takes integer k.
%   Unknown fields contain [] and are not inferred by numerical eig().
%
%   This function intentionally has no "properties" output, so Anymatrix
%   does not register it as a matrix generator.
if ~(ischar(family)||isstring(family)) || isempty(varargin)
    error('unitary_spectral_reference:Arguments', ...
          'Supply the generator name followed by its parameters.');
end
family=char(family);
n=varargin{1};
if strcmp(family,'wigner_rotation')
    validateattributes(n,{'numeric'},{'scalar','real','finite','nonnegative'});
    if 2*n~=round(2*n)
        error('unitary_spectral_reference:Spin','Wigner spin must be an integer or half-integer.');
    end
else
    validateattributes(n,{'numeric'},{'scalar','integer','finite','nonnegative'});
end
if n==0 && ~any(strcmp(family,{'jaynes_cummings','wigner_rotation'}))
    error('unitary_spectral_reference:Dimension','N must be positive.');
end
blank=struct('kind','unknown','density',[],'support',[], ...
             'atoms',[],'moments',[]);
S=struct('family',family,'dimension',double(n), ...
    'normalization','probability (density with respect to dtheta)', ...
    'finite',struct('kind','unknown','eigenvalues',[], ...
                    'moments',[],'atoms',[]), ...
    'limit',blank, ...
    'ensemble',blank, ...
    'underlying',blank, ...
    'reference','', ...
    'notes','No analytical trace spectrum is provided for this family.');

switch family
    case {'dirac_automaton','split_step_walk','driven_kitaev','julia_unitary'}
        S.dimension=2*double(n);
    case {'floquet_ising','selfdual_kicked_ising','dual_unitary_circuit','xxz_floquet'}
        S.dimension=2^double(n);
    case {'haar_unimodular_cse'}
        S.dimension=2*double(n);
    case 'entangled_reflection'
        S.dimension=double(n)^2;
    case 'wigner_rotation'
        S.dimension=2*double(n)+1;
    case 'spherical_harmonic_rotation'
        S.dimension=double(n)^2;
    case 'haar_tensor_power'
        t=param(varargin,2,2);
        S.dimension=double(n)^double(t);
    case 'haar_block_tensor_power'
        t1=param(varargin,2,1);
        t2=param(varargin,3,1);
        S.dimension=double(n)^(double(t1)+double(t2));
    case 'tensor_permutation'
        permutation=param(varargin,2,[2 1]);
        S.dimension=double(n)^numel(permutation);
    case 'kicked_bose_hubbard'
        S.dimension=double(n)+1;
    case 'jaynes_cummings'
        S.dimension=2*(double(n)+1);
    case 'chalker_coddington'
        ny=param(varargin,2,[]);
        validateattributes(ny,{'numeric'},{'scalar','integer','positive','finite'});
        S.dimension=4*double(n)*double(ny);
end

switch family
    case 'wigner_rotation'
        alpha=param(varargin,2,.3);
        beta=param(varargin,3,.7);
        gamma=param(varargin,4,.4);
        validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
        validateattributes(beta,{'numeric'},{'scalar','real','finite'});
        validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
        delta=2*acos(max(-1,min(1,cos(beta/2)*cos((alpha+gamma)/2))));
        S=exact_finite(S,@() exp(-1i*(-n:n)'*delta));
        S.reference='Gimbutas and Greengard, JCP 228 (2009)';
        S.notes='Exact SO(3) rotation eigenphases; independent of the harmonic-coordinate axis.';

    case 'spherical_harmonic_rotation'
        alpha=param(varargin,2,.3);
        beta=param(varargin,3,.7);
        gamma=param(varargin,4,.4);
        validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
        validateattributes(beta,{'numeric'},{'scalar','real','finite'});
        validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
        delta=2*acos(max(-1,min(1,cos(beta/2)*cos((alpha+gamma)/2))));
        S=exact_finite(S,@() spherical_rotation_eigs(n,delta));
        S.finite.moments=@(k) spherical_rotation_moments(n,delta,k);
        S.reference='Gimbutas and Greengard, JCP 228 (2009); derived finite character identity';
        S.notes='Exact L^2 phase multiplicities and Fejer-square trace moments. For a size-dependent delta=c/L, 0<c<pi, the limiting trace density is triangular. No limit is asserted for fixed Euler angles.';

    case 'periodic_sinc_shift'
        delta=param(varargin,2,.375);
        variant=param(varargin,3,'real');
        validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
        validateattributes(delta,{'numeric'},{'scalar','real','finite'});
        if ~(ischar(variant)||isstring(variant))
            error('unitary_spectral_reference:Variant','Invalid variant.');
        end
        variant=lower(char(variant));
        if mod(n,2)==1
            if ~strcmp(variant,'real')
                error('unitary_spectral_reference:Variant', ...
                    'Complex and Nyquist-frozen variants require even N.');
            end
            k=(-(n-1)/2:(n-1)/2)';
            S.notes='Exact centered-frequency phases of the real orthogonal odd-order fractional shift.';
        else
            k=(0:n-1)';
            k(k>n/2)=k(k>n/2)-n;
            switch variant
                case 'complex'
                    S.notes='Even-order complex unitary fractional shift with +N/2 Nyquist frequency; U(1) is the cyclic permutation.';
                case 'nyquist-frozen'
                    k(n/2+1)=0;
                    S.notes='Even-order real orthogonal fractional shift; fixed +1 Nyquist eigenvalue means U(1) is not the cyclic permutation.';
                otherwise
                    error('unitary_spectral_reference:Order', ...
                        'Even N requires an explicit valid variant.');
            end
        end
        S=exact_finite(S,@() exp(-2i*pi*k*delta/n));
        S.reference='Periodic trigonometric interpolation / normalized DFT; Puschel and Moura (2008)';

    case 'cycle_shift'
        S=exact_finite(S,@() exp(2i*pi*(0:n-1)'/n));
        S.finite.moments=@(k) shift_moments(k,n);
        S.limit=uniform_limit(S.limit);
        S.reference='Davis, Circulant Matrices (1979)';
        S.notes='Finite measure: equally weighted roots of unity; limiting trace law: uniform.';

    case 'weighted_shift'
        theta=param(varargin,2,pi*(0:n-1)'/n);
        validate_phases(theta,n);
        S=exact_finite(S,@() weighted_shift_eigs(n,theta));
        S.limit=uniform_limit(S.limit);
        S.notes='Finite eigenvalues are the roots of the product of cycle weights.';

    case {'diagonal_phase','phase_circulant'}
        supplied=numel(varargin)>=2 && ~isempty(varargin{2});
        theta=param(varargin,2,2*pi*(0:n-1)'/n);
        validate_phases(theta,n);
        S=exact_finite(S,@() exp(1i*theta(:)));
        if ~supplied
            S.limit=uniform_limit(S.limit);
            S.notes='Default equally spaced phases; uniform limiting trace law.';
        else
            S.notes='Prescribed finite phases; no limiting law follows without a size-dependent phase rule.';
        end

    case 'orthogonal_blocks'
        theta=param(varargin,2,0.3*(1:floor(n/2))');
        epsilon=param(varargin,3,1);
        validate_phases(theta,floor(n/2));
        if ~(isscalar(epsilon)&&isnumeric(epsilon)&&isreal(epsilon) ...
                && any(epsilon==[-1 1]))
            error('unitary_spectral_reference:Sign','EPSILON must be +1 or -1.');
        end
        S=exact_finite(S,@() rotation_eigs(n,theta,epsilon));
        S.notes='Exact conjugate-paired finite spectrum including odd-order real eigenvalue.';

    case 'cycle_propagator'
        t=param(varargin,2,1);
        validateattributes(t,{'numeric'},{'scalar','real','finite'});
        S=exact_finite(S,@() exp(-2i*t*cos(2*pi*(0:n-1)'/n)));
        if t==0
            S.limit.kind='atomic';
            S.limit.atoms=atoms(0,1);
        else
            a=2*abs(t);
            S.limit.kind='analytic';
            S.limit.density=@(theta) wrapped_arcsine(theta,a);
            S.limit.moments=@(k) besselj(0,2*t*k);
            if a<pi
                S.limit.support=[-a a];
            else
                S.limit.support=[-pi pi];
            end
        end
        S.reference='Ashcroft and Mermin, Solid State Physics (1976)';
        S.notes='Wrapped arcsine limiting trace law; no unfolding assumed.';

    case 'dirac_automaton'
        m=param(varargin,2,0.4);
        validateattributes(m,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
        nu=sqrt((1-m)*(1+m));
        S=exact_finite(S,@() dirac_eigs(n,nu));
        if m==0
            S.limit=uniform_limit(S.limit);
        elseif m==1
            S.limit.kind='atomic';
            S.limit.atoms=atoms([-pi/2;pi/2],[.5;.5]);
        else
            edge=asin(m);
            S.limit.kind='analytic';
            S.limit.density=@(theta) dirac_density(theta,nu);
            S.limit.support=[-pi+edge,-edge;edge,pi-edge];
        end
        S.reference='Bisio, D''Ariano and Tosini, Ann. Phys. 354 (2015)';
        S.notes='Finite trace eigenphases are exact; density has two bands and square-root edges for 0<m<1.';

    case 'grover_operator'
        marked=param(varargin,2,1);
        validateattributes(n,{'numeric'},{'scalar','integer','>=',2});
        validateattributes(marked,{'numeric'},{'scalar','integer','>=',1,'<',n});
        alpha=asin(sqrt(marked/n));
        S=exact_finite(S,@() grover_eigs(n,marked,alpha));
        S.finite.atoms=atoms([2*alpha;-2*alpha;0;pi], ...
                             [1;1;marked-1;n-marked-1]/n);
        S.reference='Grover, STOC (1996)';
        S.notes='Exact finite atomic counting measure; multiplicities are essential.';

    case 'discrete_fourier'
        S=exact_finite(S,@() dft_eigs(n));
        S.finite.atoms=dft_atoms(n);
        S.notes='F^4=I; phases 0, pi, +/-pi/2 with exact n-dependent multiplicities.';

    case 'walsh_hadamard'
        ordering=param(varargin,2,'sequency');
        if any(strcmpi(char(ordering),{'natural','hadamard'}))
            S=exact_finite(S,@() sylvester_eigs(n));
            S.finite.atoms=atoms([0;pi],[ceil(n/2);floor(n/2)]/n);
            S.notes='Natural Sylvester order is a symmetric involution.';
        else
            S.notes='Sequency row permutations generally change the transform eigenvalues.';
        end

    case 'time_frequency_shift'
        p=param(varargin,2,1);
        q=param(varargin,3,1);
        validateattributes(p,{'numeric'},{'scalar','integer','finite'});
        validateattributes(q,{'numeric'},{'scalar','integer','finite'});
        S=exact_finite(S,@() time_frequency_eigs(n,p,q));
        S.notes='Finite spectrum determined by monomial permutation cycles and phase products.';

    case 'jaynes_cummings'
        g=param(varargin,2,.8);
        delta=param(varargin,3,.35);
        t=param(varargin,4,1);
        wc=param(varargin,5,1);
        for v=[g,delta,t,wc]
            validateattributes(v,{'numeric'},{'scalar','real','finite'});
        end
        S=exact_finite(S,@() jaynes_eigs(n,g,delta,t,wc));
        S.notes='Complete truncated Hamiltonian: n excitation doublets and two boundary states.';
        S.reference='Jaynes and Cummings, Proc. IEEE 51 (1963)';

    case {'haar_special_unitary','haar_unimodular_coe','haar_unimodular_cse'}
        % Nishigaki (2024) gives exact finite-size mean spectral laws.
        % In beta=1, N=2 the special density |sin(theta)|/4 replaces
        % the generic cos(N*theta) expression; see author's footnote 1.
        S.reference='Nishigaki, PTEP 2024, 021A02, doi:10.1093/ptep/ptae018';
        S.limit=uniform_limit(S.limit);
        if n==1
            S.ensemble.kind='atomic';
            S.ensemble.atoms=atoms(0,1);
            if strcmp(family,'haar_unimodular_cse')
                S=exact_finite(S,@()ones(2,1));
            else
                S=exact_finite(S,@()1);
            end
            S.notes='The order-one unimodular ensemble is deterministic; no continuous finite-size density.';
        elseif strcmp(family,'haar_unimodular_coe') && n==2
            S.ensemble.kind='exact_expectation';
            S.ensemble.density=@(theta) abs(sin(theta))/4;
            S.ensemble.support=[-pi pi];
            S.ensemble.moments=@(k) coe_two_moments(k);
            S.notes='Exceptional unimodular COE of order two: |sin(theta)|/4 normalized density; NOT the generic beta=1 formula.';
        else
            [a,b]=unimodular_coefficients(family,n);
            S.ensemble.kind='exact_expectation';
            S.ensemble.density=@(theta) (1+a*cos(n*theta)+b*cos(2*n*theta))/(2*pi);
            S.ensemble.support=[-pi pi];
            S.ensemble.moments=@(k) unimodular_moments(k,n,a,b);
            S.notes='Exact FINITE-N ensemble mean; one random realization still has an atomic trace measure.';
        end

    case 'entangled_reflection'
        phi=param(varargin,2,pi);
        validateattributes(phi,{'numeric'},{'scalar','real','finite'});
        D=double(n)^2;
        lam=exp(1i*phi);
        if phi==pi,lam=-1;end
        if D<=2e6
            S=exact_finite(S,@() [lam;ones(D-1,1)]);
        else
            S.finite.kind='exact';
            S.finite.eigenvalues=@() entangled_eigenvalues(D,lam);
        end
        S.dimension=D;
        S.finite.moments=@(k) entangled_moments(k,D,lam);
        if abs(lam-1)<1e-14
            S.finite.atoms=atoms(0,1);
        else
            S.finite.atoms=atoms([angle(lam);0],[1/D;(D-1)/D]);
        end
        S.reference='Local certification of unitary operations, Sci. Rep. (2024)';
        S.notes='Known rank-one exceptional phase; finite atomic trace measure, not a density.';

    case 'tensor_permutation'
        perm=param(varargin,2,[2 1]);
        validateattributes(perm,{'numeric'},{'vector','integer','positive','finite'});
        perm=perm(:)';
        if ~isequal(sort(perm),1:numel(perm))
            error('unitary_spectral_reference:Permutation','Invalid tensor permutation.');
        end
        S.dimension=double(n)^numel(perm);
        S.finite.kind='exact';
        S.finite.eigenvalues=@() tensor_permutation_eigs(n,perm);
        S.finite.moments=@(k) tensor_permutation_moments(k,n,perm);
        if numel(perm)==2 && isequal(perm,[2 1])
            d=double(n);
            S.finite.atoms=atoms([0;pi],[(d+1)/(2*d);(d-1)/(2*d)]);
        end
        S.notes='Exact finite spectral moments from tensor index-cycle counts; no numerical eig() is needed.';

    case 'haar_tensor_power'
        t=param(varargin,2,2);
        seed=param(varargin,3,[]);
        validateattributes(t,{'numeric'},{'scalar','integer','positive','finite'});
        validate_tensor_seed(seed);
        V=haar_unitary(n,seed);
        z=eig(V);
        z=z./abs(z);
        S.dimension=double(n)^t;
        S.finite.kind='exact';
        S.finite.eigenvalues=@() tensor_sample_eigs(z,[],t,0);
        S.finite.moments=@(k) tensor_sample_moments(k,z,[],t,0);
        S.notes='Seeded finite Haar sample. Eigenvalues of the small base factor are computed numerically; the full tensor-product eigenvalues/moments follow exactly as products. No ensemble mean density is asserted.';

    case 'haar_block_tensor_power'
        t1=param(varargin,2,1);
        t2=param(varargin,3,1);
        seed=param(varargin,4,[]);
        validateattributes(t1,{'numeric'},{'scalar','integer','positive','finite'});
        validateattributes(t2,{'numeric'},{'scalar','integer','positive','finite'});
        validate_tensor_seed(seed);
        previous_rng=rng;
        restore=onCleanup(@()rng(previous_rng)); %#ok<NASGU>
        rng(seed,'twister');
        V1=haar_unitary(n);
        V2=haar_unitary(n);
        z1=eig(V1); z2=eig(V2);
        z1=z1./abs(z1);z2=z2./abs(z2);
        S.dimension=double(n)^(t1+t2);
        S.finite.kind='exact';
        S.finite.eigenvalues=@() tensor_sample_eigs(z1,z2,t1,t2);
        S.finite.moments=@(k) tensor_sample_moments(k,z1,z2,t1,t2);
        S.notes='Seeded independent two-base Haar sample. Only small base-factor eigenvalues require numerical eig(); full-space product reference exact.';

    case {'haar_unitary','circular_beta'}
        S.ensemble.kind='exact_expectation';
        S.ensemble.density=@(theta) ones(size(theta))/(2*pi);
        S.ensemble.support=[-pi pi];
        S.limit=uniform_limit(S.limit);
        S.notes='Uniform ensemble-mean trace law; a realized finite sample is not uniformly spaced.';

    case 'haar_orthogonal'
        S.limit=uniform_limit(S.limit);
        S.notes='Asymptotically uniform empirical phase law; finite-size O(n) parity constraints remain.';

    case 'jacobi_szego'
        a=param(varargin,2,0);
        b=param(varargin,3,0);
        mu=js_jacobi_measure(a,b);
        S.underlying.kind='analytic';
        S.underlying.density=mu.circle_density;
        S.underlying.support=[-pi pi];
        S.reference='Simon, OPUC on One Foot (2005), Sec. 11';
        S.notes='UNDERLYING is the infinite Szego-lift cyclic measure, NOT the finite CMV trace density.';

    case 'bernstein_szego'
        a=param(varargin,2,0.4);
        mu=szego_bernstein_measure(a);
        S.underlying.kind='analytic';
        S.underlying.density=mu.circle_density;
        S.underlying.support=[-pi pi];
        S.notes='Poisson/Bernstein-Szego underlying OPUC measure, not finite matrix counting density.';

    case 'circular_jacobi'
        a=param(varargin,2,0.3);
        eta=param(varargin,3,0.5);
        mu=szego_circular_jacobi_measure(a,eta);
        S.underlying.kind='analytic';
        S.underlying.density=mu.circle_density;
        S.underlying.support=[-pi pi];
        S.notes='Generalized circular Jacobi weight is the underlying OPUC vector measure.';

    case 'rogers_szego'
        q=param(varargin,2,0.6);
        mu=szego_rogers_measure(q);
        S.underlying.kind='analytic';
        S.underlying.density=mu.circle_density;
        S.underlying.support=[-pi pi];
        S.notes='Wrapped-normal underlying OPUC measure; finite CMV trace law is different.';

    case 'geronimus_cmv'
        a=param(varargin,2,0.3+0.2i);
        if ~(isnumeric(a)&&isscalar(a)&&isfinite(a)&&abs(a)<1)
            error('unitary_spectral_reference:Alpha','Coefficient must be in the open disk.');
        end
        mu=szego_geronimus_measure(a);
        S.limit.kind='support_only';
        edge=2*asin(abs(a));
        if a==0,S.limit.support=[-pi pi];
        else,S.limit.support=[-pi -edge;edge pi];end
        S.underlying.kind='analytic';
        S.underlying.density=mu.circle_density;
        S.underlying.support=mu.ac_support;
        S.underlying.atoms=mu.point_mass;
        S.underlying.poisson_density=mu.poisson_density;
        S.underlying.ac_mass=mu.ac_mass;
        S.notes=['Infinite constant-Verblunsky e1 measure has an exact ' ...
            'absolutely continuous arc density and may have an atom. ' ...
            'The finite CMV trace spectrum is atomic and its trace ' ...
            'weights 1/n are NOT the infinite e1 weights.'];

    case 'periodic_cmv'
        pattern=param(varargin,2,[.25;.2i]);
        if ~(isnumeric(pattern)&&isvector(pattern)&&~isempty(pattern) ...
                && all(isfinite(pattern(:)))&&all(abs(pattern(:))<1))
            error('unitary_spectral_reference:Alpha','Period must contain open-disk coefficients.');
        end
        if numel(pattern)==1
            edge=2*asin(abs(pattern));
            S.limit.kind='support_only';
            S.limit.support=[-pi -edge;edge pi];
        end
        S.notes='Infinite periodic CMV has band/gap structure; general band endpoints not yet encoded. Finite boundary eigenvalues may lie in gaps.';

    otherwise
        % No analytic trace information is asserted for this family.
end
end

function v=param(args,k,default)
if numel(args)<k || isempty(args{k})
    v=default;
else
    v=args{k};
end
end

function S=exact_finite(S,fn)
S.finite.kind='exact';
S.finite.eigenvalues=fn;
S.finite.moments=@(k) finite_moments(fn,k);
end

function m=finite_moments(fn,k)
validateattributes(k,{'numeric'},{'integer','finite','real'});
z=fn();
sz=size(k);
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),sz);
end

function m=shift_moments(k,n)
validateattributes(k,{'numeric'},{'integer','finite','real'});
m=double(mod(k,n)==0);
end

function d=uniform_limit(d)
d.kind='analytic';
d.density=@(theta) ones(size(theta))/(2*pi);
d.support=[-pi pi];
d.moments=@(k) double(k==0);
end

function A=atoms(theta,weights)
A=struct('angles',theta(:),'weights',weights(:));
end

function validate_phases(theta,n)
if ~(isnumeric(theta)&&isvector(theta)&&isreal(theta) ...
        && all(isfinite(theta(:)))&&numel(theta)==n)
    error('unitary_spectral_reference:Phases', ...
          'THETA must be a finite real vector with the required number of entries.');
end
end

function z=weighted_shift_eigs(n,theta)
phi=angle(exp(1i*sum(theta(:))));
z=exp(1i*(phi+2*pi*(0:n-1)')/n);
end

function z=rotation_eigs(n,theta,epsilon)
z=reshape([exp(1i*theta(:)).';exp(-1i*theta(:)).'],[],1);
if mod(n,2),z=[z;epsilon];end
end

function z=dirac_eigs(n,nu)
omega=acos(nu*cos(2*pi*(0:n-1)'/n));
z=[exp(1i*omega);exp(-1i*omega)];
end

function rho=dirac_density(theta,nu)
t=angle(exp(1i*theta));
disc=nu^2-cos(t).^2;
rho=zeros(size(t));
interior=disc>0;
rho(interior)=abs(sin(t(interior)))./(2*pi*sqrt(disc(interior)));
edge=abs(disc)<=16*eps(max(1,nu^2));
rho(edge & ~interior)=Inf;
end

function rho=wrapped_arcsine(theta,a)
t=angle(exp(1i*theta));
rho=zeros(size(t));
K=ceil((a+pi)/(2*pi));
for j=-K:K
    x=t+2*pi*j;
    d=a^2-x.^2;
    inside=d>0;
    rho(inside)=rho(inside)+1./(pi*sqrt(d(inside)));
    onedge=abs(abs(x)-a)<=16*eps(max(1,a));
    rho(onedge)=Inf;
end
end

function z=grover_eigs(n,marked,alpha)
z=[exp(2i*alpha);exp(-2i*alpha); ...
   ones(marked-1,1);-ones(n-marked-1,1)];
end

function counts=dft_counts(n)
switch mod(n,4)
    case 0, t=1-1i;
    case 1, t=1;
    case 2, t=0;
    otherwise, t=-1i;
end
fixed=1+double(mod(n,2)==0);
counts=round([(n+fixed+2*real(t))/4; ...
              (n+fixed-2*real(t))/4; ...
              (n-fixed+2*imag(t))/4; ...
              (n-fixed-2*imag(t))/4]);
end

function z=dft_eigs(n)
z=repelem([1;-1;1i;-1i],dft_counts(n));
z=z(:);
end

function a=dft_atoms(n)
a=atoms([0;pi;pi/2;-pi/2],dft_counts(n)/n);
end

function z=sylvester_eigs(n)
z=[ones(ceil(n/2),1);-ones(floor(n/2),1)];
end

function z=time_frequency_eigs(n,p,q)
rows=mod((0:n-1)+q,n)+1;
phase=exp(2*pi*1i*mod(p,n)*(rows-1)/n);
seen=false(1,n);
z=complex(zeros(n,1));used=0;
for j=1:n
    if seen(j),continue;end
    idx=j; phi=0;len=0;
    while ~seen(idx)
        seen(idx)=true;
        len=len+1;
        phi=phi+angle(phase(idx));
        idx=rows(idx);
    end
    z(used+(1:len))=exp(1i*(phi+2*pi*(0:len-1)')/len);
    used=used+len;
end
end

function z=jaynes_eigs(n,g,delta,t,wc)
wa=wc+delta;
k=(0:n-1)';
center=wc*(k+.5);
splitting=hypot(delta/2,g*sqrt(k+1));
energy=[-wa/2;center+splitting;center-splitting;wc*n+wa/2];
z=exp(-1i*t*energy);
end


function [a,b]=unimodular_coefficients(family,n)
% Fourier coefficients relative to 1/(2*pi) in Nishigaki theorem (1).
a=0;b=0;
switch family
    case 'haar_special_unitary'
        a=-(-1)^n*(2/n);
    case 'haar_unimodular_coe'
        loga=0.5*log(pi)+gammaln(n)-(n-1)*log(2) ...
             -gammaln(n/2+1.5)-gammaln(n/2+1);
        a=-(-1)^n*exp(loga);
    case 'haar_unimodular_cse'
        loga=2*n*log(2)+2*gammaln(n+1) ...
             -gammaln(2*n+1)-log(n);
        a=-(-1)^n*exp(loga);
        b=2/((2*n-1)*n);
end
end

function m=unimodular_moments(k,n,a,b)
validateattributes(k,{'numeric'},{'integer','finite','real'});
m=double(k==0)+(a/2)*double(abs(k)==n)+(b/2)*double(abs(k)==2*n);
end

function m=coe_two_moments(k)
validateattributes(k,{'numeric'},{'integer','finite','real'});
m=zeros(size(k));
m(k==0)=1;
mask=(k~=0 & mod(k,2)==0);
m(mask)=-1./(k(mask).^2-1);
end

function validate_tensor_seed(seed)
if isempty(seed)
    error('unitary_spectral_reference:Seed', ...
        'Seeded Haar tensor samples require an explicit nonnegative SEED.');
end
validateattributes(seed,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
end

function z=entangled_eigenvalues(D,lambda)
if D>2e6
    error('unitary_spectral_reference:Materialization', ...
          'Exact eigenvalue enumeration limited to 2e6; use finite.moments.');
end
z=[lambda;ones(D-1,1)];
end

function m=entangled_moments(k,D,lambda)
validateattributes(k,{'numeric'},{'integer','finite','real'});
m=((D-1)+lambda.^k)/D;
end

function m=tensor_permutation_moments(k,d,perm)
validateattributes(k,{'numeric'},{'integer','finite','real'});
nFactors=numel(perm);
cycleLengths=permutation_cycle_lengths(perm);
m=zeros(size(k));
for j=1:numel(k)
    ncycles=sum(gcd(cycleLengths,abs(k(j))));
    m(j)=d^(ncycles-nFactors);
end
end

function lens=permutation_cycle_lengths(perm)
t=numel(perm);
seen=false(1,t);
lens=[];
for i=1:t
    if seen(i),continue;end
    j=i;ell=0;
    while ~seen(j)
        seen(j)=true;
        ell=ell+1;
        j=perm(j);
    end
    lens(end+1)=ell; %#ok<AGROW>
end
end

function z=tensor_permutation_eigs(d,perm)
row=um_tensor_permutation_map(d,perm);
n=numel(row);
seen=false(n,1);
z=complex(zeros(n,1));
used=0;
for i=1:n
    if seen(i),continue;end
    j=i;ell=0;
    while ~seen(j)
        seen(j)=true;
        ell=ell+1;
        j=row(j);
    end
    z(used+(1:ell))=exp(2i*pi*(0:ell-1)'/ell);
    used=used+ell;
end
end

function z=tensor_sample_eigs(z1,z2,t1,t2)
n=numel(z1)^t1;
if t2>0,n=n*numel(z2)^t2;end
if n>2e6
    error('unitary_spectral_reference:Materialization', ...
        'Exact eigenvalue enumeration limited to 2e6; use finite.moments.');
end
z=1;
for j=1:t1,z=kron(z,z1);end
for j=1:t2,z=kron(z,z2);end
z=z(:);
end

function m=tensor_sample_moments(k,z1,z2,t1,t2)
validateattributes(k,{'numeric'},{'integer','finite','real'});
shape=size(k);
kk=reshape(k,1,[]);
m=mean(bsxfun(@power,z1(:),kk),1).^t1;
if t2>0
    m=m.*mean(bsxfun(@power,z2(:),kk),1).^t2;
end
m=reshape(m,shape);
end

function z=spherical_rotation_eigs(L,delta)
m=(-(L-1):(L-1))';
z=repelem(exp(-1i*m*delta),L-abs(m));
z=z(:);
end

function M=spherical_rotation_moments(L,delta,k)
validateattributes(k,{'numeric'},{'integer','finite','real'});
shape=size(k);
kk=k(:);
M=zeros(size(kk));
j=(0:L-1)';
for r=1:numel(kk)
    M(r)=abs(sum(exp(-1i*kk(r)*delta*j)))^2/L^2;
end
M=reshape(M,shape);
end
