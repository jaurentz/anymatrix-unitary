function [U,properties] = agassi_propagator(j,epsilon,g,V,t,particles,mode)
%AGASSI_PROPAGATOR Number-conserving nuclear pairing/particle-hole model.
%   U=AGASSI_PROPAGATOR(J,EPSILON,G,V,T) constructs exp(-1i*T*H)
%   in the full Fock basis of 4J fermion orbitals (two single-particle
%   levels, each with 2J magnetic substates). J is a positive integer;
%   the four-orbital (J=1) case is the published 16-state simulation.
%
%   The Agassi Hamiltonian (Perez-Fernandez et al., Phys. Lett. B 829,
%   137133 (2022), Eqs. 1-6) is
%       H=EPSILON*J0-G*(Ap+Am)*(Ap+Am)' ...
%              -V/2*(Jp^2+(Jp')^2),
%   Jp=sum_m c_(+,m)'*c_(-,m), J0=(Nplus-Nminus)/2,
%   Ap=sum_(m>0) c_(+,m)'*c_(+,-m)' and
%   Am=sum_(m>0) c_(-,m)'*c_(-,-m)'.
%   The ordering is upper +m, upper -m, lower +m, lower -m.
%   Fermionic signs use the occupation-number Jordan-Wigner basis.
%   Different Pauli occupancy conventions may shift the energy origin
%   or change matrix entries by a fixed basis conjugation. The physical
%   eigenvalue differences and dynamics in a fixed number sector agree.
%
%   An optional PARTICLES integer restricts to the invariant sector
%   with exactly that many fermions, yielding order nchoosek(4*J,
%   PARTICLES). The default [] uses all 2^(4J) occupation states.
%   'hamiltonian' returns the sparse Hermitian H, not a unitary.
%   'spectral' gives exact finite eigenphases if G=V=0 or T=0.
%   Otherwise finite spectral data is left unknown. The unitary is
%   explicitly assembled for J<=2 (up to 256 states); 'hamiltonian'
%   and 'spectral' also support J=3.
%
%   Nuclear pairing and particle-hole physics differs from the
%   collective two-boson extended Lipkin model and Floquet spin chains.
%   https://doi.org/10.1016/j.physletb.2022.137133
properties={'square','unitary','normal','complex','parameter-dependent', ...
            'scalable'};
if nargin==0,U=[];return;end
validateattributes(j,{'numeric'},{'scalar','integer','positive','finite','<=',3});
if nargin<2 || isempty(epsilon),epsilon=1;end
if nargin<3 || isempty(g),g=.3;end
if nargin<4 || isempty(V),V=.4;end
if nargin<5 || isempty(t),t=1;end
for v=[epsilon g V t]
    validateattributes(v,{'numeric'},{'scalar','real','finite'});
end
if nargin<6,particles=[];end
if ~isempty(particles)
    validateattributes(particles,{'numeric'}, ...
        {'scalar','integer','finite','>=',0,'<=',4*j});
end
if nargin<7 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('agassi_propagator:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','hamiltonian','spectral'}))
    error('agassi_propagator:Mode', ...
        'Supported modes are matrix, hamiltonian, spectral.');
end
if strcmp(mode,'matrix') && j>2
    error('agassi_propagator:ExplicitOrder', ...
        'For J=3 use hamiltonian mode and a structured exponential method.');
end
nf=4*j;
allstates=uint32((0:2^nf-1)');
if isempty(particles)
    take=(1:numel(allstates))';
else
    count=zeros(numel(allstates),1);
    for p=1:nf
        count=count+double(bitget(allstates,p));
    end
    take=find(count==particles);
end
N=numel(take);
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',N);
    U.family='agassi_propagator';
    U.reference='Perez-Fernandez et al., Phys. Lett. B 829 (2022), 137133';
    U.notes=['Agassi pairing and particle-hole fermionic dynamics. ' ...
       'An optional number-conserving sector reduces Fock dimension. ' ...
       'Generic finite eigenphase density is not known analytically.'];
    if t==0
        zfun=@()ones(N,1);
    elseif g==0 && V==0
        en=zeros(N,1);
        sel=allstates(take);
        for p=1:2*j
            en=en+.5*epsilon*double(bitget(sel,p));
        end
        for p=2*j+1:4*j
            en=en-.5*epsilon*double(bitget(sel,p));
        end
        zfun=@()exp(-1i*t*en);
    else
        zfun=[];
    end
    if ~isempty(zfun)
        U.finite.kind='exact';
        U.finite.eigenvalues=zfun;
        U.finite.moments=@(k)ap_finite_moments(zfun,k);
    end
    return
end
d=double(numel(allstates));
C=cell(nf,1);
for p=1:nf
    occupied=find(bitget(allstates,p)~=0);
    rows=occupied-2^(p-1);
    if p==1
        parity=zeros(numel(occupied),1);
    else
        parity=zeros(numel(occupied),1);
        for q=1:p-1
            parity=parity+double(bitget(allstates(occupied),q));
        end
    end
    values=1-2*mod(parity,2);
    C{p}=sparse(rows,occupied,values,d,d);
end
Jp=sparse(d,d);J0=sparse(d,d);
Ap=sparse(d,d);Am=sparse(d,d);
for p=1:2*j
    Jp=Jp+C{p}'*C{2*j+p};
    J0=J0+.5*(C{p}'*C{p}-C{2*j+p}'*C{2*j+p});
end
for p=1:j
    Ap=Ap+C{p}'*C{j+p}';
    Am=Am+C{2*j+p}'*C{3*j+p}';
end
Pair=Ap+Am;
H=epsilon*J0-g*(Pair*Pair')-.5*V*(Jp*Jp+Jp'*Jp');
H=H(take,take);
if strcmp(mode,'hamiltonian')
    U=H;
    return
end
U=expm(-1i*t*full(H));
end

function m=ap_finite_moments(zfun,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
z=zfun();
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),size(k));
end
