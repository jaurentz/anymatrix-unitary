function [U,properties]=bosonic_symmetric_power(d,p,V,output)
%BOSONIC_SYMMETRIC_POWER Number-conserving bosonic second quantization.
%   U=BOSONIC_SYMMETRIC_POWER(D,P,V) lifts a D-by-D unitary V to the
%   P-boson sector Sym^P(C^D), order N=nchoosek(D+P-1,P).
%   Basis occupation vectors (n_1,...,n_D) sum to P, including repeated
%   occupation. The ordering is lexicographic in bars-and-stars separators.
%   U(m,n)=per(V[m repeated rows, n repeated columns])/
%          sqrt(prod(factorial(m_j))*prod(factorial(n_j))).
%   This is strictly unitary by the bosonic Fock-space representation.
%   Default V=normalized DFT(D), P=min(2,D). P=0 gives [1],
%   P=1 gives V. Small explicit matrices only (N<=100, P<=7).
%   OUTPUT='spectral' gives EXACT N eigenphases as all occupation
%   products of eig(V), counting multiplicities, without forming U.
%   Fourier moments are complete homogeneous symmetric polynomials
%   h_P(eig(V).^k)/N, computed in O(D*P) work per k.
%   This does NOT claim a fast matrix-free general permanent action,
%   and is not the full D^P tensor power nor the fermionic wedge^P.
%   [117] S. Scheel (2004) https://arxiv.org/abs/quant-ph/0406127
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(p),p=min(2,d);end
validateattributes(p,{'numeric'},{'scalar','integer','nonnegative','finite'});
if nargin<3||isempty(V),V=discrete_fourier(d);end
if ~isnumeric(V)||~isequal(size(V),[d d])||any(~isfinite(V(:)))
    error('bosonic_symmetric_power:Base','V must be finite D-by-D.');
end
if norm(V'*V-eye(d),'fro')>2e-10*max(1,d)
    error('bosonic_symmetric_power:Unitary','V must be unitary.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('bosonic_symmetric_power:Mode','OUTPUT must be scalar text.');
end
N=nchoosek(d+p-1,p);
if ~isfinite(N)||N>flintmax
    error('bosonic_symmetric_power:Order','Dimension exceeds exact index range.');
end
switch lower(char(output))
    case 'spectral'
        z=eig(full(V));
        z=z./abs(z);
        U=unitary_spectral_reference('bosonic_symmetric_power',d);
        U.dimension=N;
        U.finite.kind='exact';
        U.finite.eigenvalues=@()um_boson_spectral(z,p,'eigenvalues');
        U.finite.moments=@(k)um_boson_spectral(z,p,'moments',k);
        U.reference='Scheel (2004), Permanents in linear optical networks';
        U.notes=['Normalized P-boson Fock trace phases are all repeated ' ...
            'occupation products of one-particle eigenvalues.'];
    case {'matrix','explicit'}
        if p==0,U=1;return;end
        if p==1,U=V;return;end
        if d==1,U=V(1,1)^p;return;end
        if N>100||p>7
            error('bosonic_symmetric_power:Order', ...
                'Explicit permanents require N<=100 and P<=7; use spectral.');
        end
        occ=um_boson_occupations(d,p);
        U=zeros(N,N);
        if ~isreal(V),U=complex(U);end
        repeats=cell(N,1);
        normfactor=ones(N,1);
        for a=1:N
            repeats{a}=repelem(1:d,occ(a,:));
            normfactor(a)=sqrt(prod(factorial(occ(a,:))));
        end
        for a=1:N
            for b=1:N
                sub=V(repeats{a},repeats{b});
                U(a,b)=um_boson_permanent(sub)/(normfactor(a)*normfactor(b));
            end
        end
    otherwise
        error('bosonic_symmetric_power:Mode', ...
            'Only matrix and spectral modes are supported.');
end
end
