function S = ep_spectral_reference(family,n,varargin)
%EP_SPECTRAL_REFERENCE Standard analytical trace descriptors for pencil lifts.
%   Does not construct a matrix or advance the random number generator.
%   S.FINITE is the equal-weight trace measure (never OPUC vector weights).
switch family
    case 'mobius_unitary_pencil'
        a=varargin{1};
        S=unitary_spectral_reference('cycle_shift',n);
        zfun=@() (exp(2i*pi*(0:n-1)'/n)-a) ./ ...
            (1-conj(a)*exp(2i*pi*(0:n-1)'/n));
        S.finite.kind='exact';
        S.finite.eigenvalues=zfun;
        S.finite.moments=@(k) ep_finite_moments(zfun,k);
        S.limit.kind='analytic';
        S.limit.density=@(theta) (1-abs(a)^2) ./ ...
            (2*pi*abs(1+conj(a)*exp(1i*theta)).^2);
        S.limit.support=[-pi pi];
        S.limit.moments=@(k) ep_poisson_moments(k,a);
        S.notes='Exact disk-Mobius images of cyclic roots; Poisson limiting trace density.';
        S.reference='Simon, Orthogonal Polynomials on the Unit Circle (2005)';
    case 'left_weighted_unitary_pencil'
        S=unitary_spectral_reference('cycle_shift',n);
        S.notes='Exactly the cyclic-shift spectrum despite ill-scaled pencil coefficients.';
    case 'unitary_delay_pencil'
        delays=varargin{1};
        S=unknown_descriptor(n,sum(delays));
        S.notes='Integer-delay network unitary lift; roots of the nonlinear secular determinant.';
        S.reference='Kottos and Smilansky, Ann. Phys. 274 (1999)';
    case 'block_cmv_pencil'
        q=varargin{1};
        S=unknown_descriptor(n,n*q);
        S.notes='Unitary block-CMV closure; finite trace law depends on block coefficients.';
        S.reference='Damanik, Pushnitski and Simon, Surveys Approx. Theory 4 (2008)';
    case 'floquet_factor_pencil'
        theta1=varargin{1};theta2=varargin{2};
        S=unitary_spectral_reference('split_step_walk',n,theta1,theta2);
        S.notes='Same spectrum as the split-step walk; generalized factor pencil.';
    case 'cayley_hermitian_pencil'
        a=varargin{1};b=varargin{2};t=varargin{3};
        S=unitary_spectral_reference('cayley_jacobi',n,a,b,t);
        S.notes='Cayley image of the Jacobi recurrence; explicit generalized pencil.';
    otherwise
        error('ep_spectral_reference:Family','Unknown pencil family %s.',family);
end
S.family=family;
end

function S=unknown_descriptor(n,N)
S=unitary_spectral_reference('cycle_shift',max(1,n));
S.dimension=double(N);
S.finite=struct('kind','unknown','eigenvalues',[], ...
    'moments',[],'atoms',[]);
S.limit=struct('kind','unknown','density',[], ...
    'support',[],'atoms',[],'moments',[]);
S.ensemble=S.limit;
S.underlying=S.limit;
end

function m=ep_finite_moments(zfun,k)
validateattributes(k,{'numeric'},{'integer','real','finite'});
z=zfun(); sz=size(k);
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),sz);
end

function m=ep_poisson_moments(k,a)
validateattributes(k,{'numeric'},{'integer','real','finite'});
m=zeros(size(k));p=k>=0;
m(p)=(-a).^k(p);
m(~p)=conj((-a).^(-k(~p)));
end
