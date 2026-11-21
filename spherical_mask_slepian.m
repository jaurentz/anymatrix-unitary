function [U,properties]=spherical_mask_slepian(L,components,mode)
%SPHERICAL_MASK_SLEPIAN Exact finite bandlimited weighted spherical masks.
%   U=SPHERICAL_MASK_SLEPIAN(L,COMPONENTS) is a real orthogonal N-by-N
%   ANALYSIS matrix Q', N=L^2, where columns Q diagonalize the real
%   spherical harmonic concentration matrix of a nonnegative mask
%      W(Omega)=sum_j w_j*1_{R_j A_j}(Omega), 0<=W<=1.
%   A_j is a north polar cap OR a full-colatitude spherical lune.
%   R_j is a proper SO(3) z-y-z Euler rotation. This is a finite
%   *weighted mixture of rotated exact building blocks*, NOT generic
%   polygon quadrature nor a claim of an exact binary union when the
%   component regions overlap.
%
%   COMPONENTS is a struct array, with equal fields:
%       type       'cap' or 'lune'
%       weight     real in [0,1], sum(weights)<=1 (overlap-safe)
%       angles     finite real [alpha beta gamma] in radians
%       parameters scalar CAP in [0,pi], or [PHISTART,WIDTH] in
%                  radians (WIDTH in [0,2*pi]).
%   For default input [] use an empty zero mask (U=eye(N)).
%   A component with nonzero angles is actively rotated in space.
%
%   The concentration in the shared real harmonic basis of
%   SPHERICAL_SLEPIAN / SPHERICAL_LUNE_SLEPIAN is
%       K=sum_j w_j Dreal(R_j)*K_base_j*Dreal(R_j)'.
%   Each K_base is exactly integrated in mathematical real
%   arithmetic: Gaussian radial quadrature and analytic longitude
%   factors. Dreal is the exact degree-block SO(3) representation
%   in real harmonic coordinates. The basis permutation and the
%   Condon--Shortley sine/cosine signs are explicitly corrected
%   from the degree-major UM_WIGNER_REAL_BASIS ordering.
%   Thus K is real symmetric PSD with 0<=K<=I in exact arithmetic.
%   Concentration factors lambda are in [0,1], but these are NOT
%   the eigenphases of orthogonal U. General eigenphases UNKNOWN.
%   Exact trace K = L^2*sum_j weight_j*area(A_j)/(4*pi).
%
%   'matrix'/'analysis': Q' (default)
%   'basis': Q
%   'concentration': nonunitary symmetric K
%   'parameters': K,Q,lambda,mask components, Shannon number, labels
%   'spectral': unknown finite eigenphases of Q', NOT lambda.
%   Empty mask uses Q=I. Nearly repeated eigenspaces are
%   canonically oriented using projected coordinate pivots and
%   a documented finite eigengap tolerance.
%
%   Explicit dense reference only, 1<=L<=12, <=8 components.
%   Rotated components cost O(L^4) storage for their full
%   representation and O(L^6) dense algebra each. The final eig(K)
%   also costs O(L^6). No subquadratic action claimed.
%
%   [156] Simons, Dahlen, Wieczorek, SIAM Rev. 48 (2006), 504--536.
%   https://doi.org/10.1137/S0036144504445765
properties={'square','unitary','orthogonal','real','normal', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','real','integer','positive','finite','<=',12});
if nargin<2||isempty(components),components=struct( ...
        'type',{},'weight',{},'angles',{},'parameters',{});end
if nargin<3||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('spherical_mask_slepian:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','analysis','basis', ...
        'concentration','parameters','spectral'}))
    error('spherical_mask_slepian:Mode','Unknown output mode.');
end
if ~isstruct(components)||numel(components)>8
    error('spherical_mask_slepian:Components', ...
        'COMPONENTS must be a struct array with at most eight elements.');
end
n=L^2;
weightSum=0;areaFraction=0;
for j=1:numel(components)
    c=components(j);
    if ~all(isfield(c,{'type','weight','angles','parameters'}))
        error('spherical_mask_slepian:Components', ...
            'Each component needs type, weight, angles, parameters.');
    end
    if ~(ischar(c.type)||(isstring(c.type)&&isscalar(c.type)))
        error('spherical_mask_slepian:Type', ...
            'Type must be cap or lune.');
    end
    if ~(isnumeric(c.weight)&&isscalar(c.weight)&& ...
            isreal(c.weight)&&isfinite(c.weight) ...
            &&c.weight>=0&&c.weight<=1)
        error('spherical_mask_slepian:Weight', ...
            'Each weight must be between zero and one.');
    end
    if ~(isnumeric(c.angles)&&isreal(c.angles) ...
            &&isvector(c.angles)&&numel(c.angles)==3 ...
            &&all(isfinite(c.angles(:))))
        error('spherical_mask_slepian:Angles', ...
            'Rotation angles must be a finite three-vector.');
    end
    switch lower(char(c.type))
        case 'cap'
            if ~(isnumeric(c.parameters)&&isscalar(c.parameters) ...
                    &&isreal(c.parameters)&&isfinite(c.parameters) ...
                    &&c.parameters>=0&&c.parameters<=pi)
                error('spherical_mask_slepian:Region', ...
                    'Cap colatitude must be in [0,pi].');
            end
            fraction=(1-cos(c.parameters))/2;
        case 'lune'
            p=c.parameters;
            if ~(isnumeric(p)&&isreal(p)&&isvector(p) ...
                    &&numel(p)==2&&all(isfinite(p(:))) ...
                    &&p(2)>=0&&p(2)<=2*pi&&abs(p(1))<=1e6)
                error('spherical_mask_slepian:Region', ...
                    'Lune needs [start,width], width in [0,2*pi].');
            end
            fraction=p(2)/(2*pi);
        otherwise
            error('spherical_mask_slepian:Type', ...
                'Only cap or lune regions are supported.');
    end
    weightSum=weightSum+c.weight;
    areaFraction=areaFraction+c.weight*fraction;
end
if weightSum>1+1e-13
    error('spherical_mask_slepian:WeightSum', ...
        'Sum of nonnegative weights must not exceed one.');
end
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('spherical_mask_slepian',n);
    U.bandlimit=L;U.components=components;
    U.shannon_number=n*areaFraction;
    U.reference='Simons--Dahlen--Wieczorek, SIAM Review 48 (2006)';
    U.notes=['Finite real orthogonal Slepian basis for exact weighted ', ...
        'mixture of rotated caps/lunes. Matrix eigenphases ', ...
        'UNKNOWN; concentration eigenvalues of nonunitary K ', ...
        'are different spectral data.'];
    return
end
% Translation between (m,cos/sin,ell) and degree-major (ell,m_real)
% harmonic orders. The real Wigner basis uses R_cos=(-1)^m
% times our Condon--Shortley cosine, R_sin=(-1)^(m+1)
% times our Condon--Shortley sine.
lookup=zeros(n,1);signs=zeros(n,1);
degree=zeros(n,1);order=zeros(n,1);kind=zeros(n,1);
h=0;
for m=0:L-1
    if m==0,types=0;else,types=[1 2];end
    for kindIndex=types
        for ell=m:L-1
            h=h+1;
            degree(h)=ell;order(h)=m;kind(h)=kindIndex;
            if kindIndex==0,mm=0;signs(h)=1;
            elseif kindIndex==1,mm=m;signs(h)=(-1)^m;
            else,mm=-m;signs(h)=(-1)^(m+1);end
            lookup(h)=ell^2+ell+mm+1;
        end
    end
end
if h~=n||numel(unique(lookup))~=n
    error('spherical_mask_slepian:Internal', ...
        'Real harmonic basis matching failed.');
end
K=zeros(n);
for j=1:numel(components)
    c=components(j);
    if c.weight==0,continue;end
    switch lower(char(c.type))
        case 'cap'
            Kbase=spherical_slepian(L,c.parameters,'concentration');
        case 'lune'
            Kbase=spherical_lune_slepian(L,c.parameters(1), ...
                c.parameters(2),'concentration');
    end
    a=c.angles(:);
    if all(a==0)
        R=eye(n);
    else
        D=full(real_spherical_harmonic_rotation( ...
            L,a(1),a(2),a(3)));
        R=(signs*signs').*D(lookup,lookup);
    end
    K=K+c.weight*(R*Kbase*R');
end
K=(K+K')/2;
if isempty(components)||weightSum==0
    Q=eye(n);lambda=zeros(n,1);perm=(1:n)';
    tieTolerance=0;
else
    [V,D]=eig(K);
    [lambda,perm]=sort(real(diag(D)),'descend');
    Q=real(V(:,perm));
    tieTolerance=2e-12*max(1,norm(K,2));
    i=1;
    while i<=n
        j=i;
        while j<n&&abs(lambda(j)-lambda(j+1))<=tieTolerance
            j=j+1;
        end
        if j>i
            P=Q(:,i:j)*Q(:,i:j)';
            W=zeros(n,j-i+1);count=0;
            for pivot=1:n
                v=P(:,pivot);
                if count>0
                    v=v-W(:,1:count)*(W(:,1:count)'*v);
                    v=v-W(:,1:count)*(W(:,1:count)'*v);
                end
                nv=norm(v);
                if nv>1e-8
                    count=count+1;W(:,count)=v/nv;
                    if count==j-i+1,break;end
                end
            end
            if count~=j-i+1
                error('spherical_mask_slepian:Degeneracy', ...
                    'Repeated-subspace canonicalization failed.');
            end
            Q(:,i:j)=W;
        end
        i=j+1;
    end
    for k=1:n
        [~,pivot]=max(abs(Q(:,k)));
        if Q(pivot,k)<0,Q(:,k)=-Q(:,k);end
    end
end
switch mode
    case {'matrix','analysis'}
        U=Q.';
    case 'basis'
        U=Q;
    case 'concentration'
        U=K;
    case 'parameters'
        U=struct('concentration',K,'basis',Q,'analysis',Q.', ...
            'concentration_eigenvalues',lambda, ...
            'bandlimit',L,'components',components, ...
            'shannon_number',n*areaFraction, ...
            'weight_sum',weightSum,'basis_permutation',lookup, ...
            'basis_signs',signs,'degree',degree, ...
            'azimuthal_order',order,'harmonic_kind',kind, ...
            'eigenvalue_original_index',perm(:), ...
            'eigenvalue_tie_tolerance',tieTolerance, ...
            'notes',['Exact in real arithmetic for weighted sums ', ...
                'of rotated caps/lunes only; overlapping regions ', ...
                'give a weighted mask, not a binary union.']);
end
end
