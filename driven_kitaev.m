function [U,properties] = driven_kitaev(n,opts,output)
%DRIVEN_KITAEV Two-step driven spinless p-wave Kitaev chain (BdG).
%   U = exp(-1i*tau2*H_BdG(mu2))*exp(-1i*tau1*H_BdG(mu1)).
%   n is number of lattice sites; U is 2*n-by-2*n in the Nambu
%   basis (c_1,...,c_n,c_1^dagger,...,c_n^dagger).
%   It is NOT the 2^n-by-2^n many-body Fock-space propagator.
%   Options struct fields: mu1 (default 0.3), mu2 (2.7),
%   hopping (1), pairing (1), tau1/tau2 (0.7), boundary ('open').
%   Real hopping, pairing and chemical potentials give BdG
%   particle-hole symmetry: tau_x*conj(U)*tau_x=U.
%   Open boundaries can support topological edge quasienergies; they
%   are not guaranteed for every choice of drive parameters.
%   References: D. J. Yates and A. Mitra, Phys. Rev. B 96 (2017),
%   115108; M. van Caspel et al., SciPost Phys. 6 (2019), 026.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [33] https://doi.org/10.1103/PhysRevB.96.115108
%   [34] https://doi.org/10.21468/SciPostPhys.6.2.026
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3
    if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
        U=unitary_spectral_reference('driven_kitaev',n,opts);
        return
    end
    error('driven_kitaev:Mode','Only terminal spectral mode is supported for this generator.');
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(opts),opts=struct();end
if ~isstruct(opts) || numel(opts)~=1
    error('driven_kitaev:Options','Options must be a scalar struct.');
end
defaults=struct('mu1',0.3,'mu2',2.7,'hopping',1,'pairing',1, ...
                'tau1',0.7,'tau2',0.7,'boundary','open');
names=fieldnames(opts);
for k=1:numel(names)
    if ~isfield(defaults,names{k})
        error('driven_kitaev:Options','Unknown option: %s.',names{k});
    end
    defaults.(names{k})=opts.(names{k});
end
p=defaults;
for k={'mu1','mu2','hopping','pairing','tau1','tau2'}
    validateattributes(p.(k{1}),{'numeric'},{'scalar','real','finite'});
end
if p.tau1<0 || p.tau2<0
    error('driven_kitaev:Durations','Step durations must be nonnegative.');
end
p.boundary=lower(char(p.boundary));
if ~any(strcmp(p.boundary,{'open','periodic'}))
    error('driven_kitaev:Boundary','Boundary must be open or periodic.');
end
if n==2 && strcmp(p.boundary,'periodic')
    error('driven_kitaev:Boundary','Periodic chains require n>=3.');
end
H1=local_bdg(n,p.mu1,p.hopping,p.pairing,p.boundary);
H2=local_bdg(n,p.mu2,p.hopping,p.pairing,p.boundary);
U=expm(-1i*p.tau2*H2)*expm(-1i*p.tau1*H1);
end

function H=local_bdg(n,mu,t,delta,boundary)
h=-mu*eye(n);
B=zeros(n);
for j=1:n-1
    h(j,j+1)=-t;
    h(j+1,j)=-t;
    B(j,j+1)=delta;
    B(j+1,j)=-delta;
end
if strcmp(boundary,'periodic')
    h(n,1)=-t;
    h(1,n)=-t;
    B(n,1)=delta;
    B(1,n)=-delta;
end
H=[h,B;B',-h.'];
end
