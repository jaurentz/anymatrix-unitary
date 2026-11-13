function S = unitary_spectral_reference(family,varargin)
%UNITARY_SPECTRAL_REFERENCE Analytic spectral reference data for unitary generators.
%   S=UNITARY_SPECTRAL_REFERENCE(FAMILY,PARAMETERS...) returns a standardized
%   descriptor, without constructing a matrix or consuming random numbers.
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
validateattributes(n,{'numeric'},{'scalar','integer','finite','nonnegative'});
if n==0 && ~strcmp(family,'jaynes_cummings')
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
            S.limit.moments=@(k) (-1i).^k.*besselj(k,2*t);
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
        S.limit.kind='support_only';
        edge=2*asin(abs(a));
        S.limit.support=[-pi -edge; edge pi];
        S.notes='Essential bands for constant Verblunsky coefficient; finite boundary phases may place eigenvalues in gaps. Underlying OPUC measure may have an atom.';

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
