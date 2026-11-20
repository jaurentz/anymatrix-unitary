function [U,properties]=fermionic_exterior_power(d,p,V,output)
%FERMIONIC_EXTERIOR_POWER Number-conserving fermionic unitary lift.
%   U=FERMIONIC_EXTERIOR_POWER(D,P,V) is wedge^P(V) for a D-by-D
%   one-particle unitary V. In increasing occupied-mode subset order,
%   U(I,J)=det(V(I,J)); 0<=P<=D, N=nchoosek(D,P).
%   Thus U is N-by-N unitary. P=0 gives [1]; P=1 gives V; P=D
%   gives [det(V)]. Default V=normalized DFT(D), default P=min(2,D).
%   'spectral' gives the exact P-mode subset products of eig(V),
%   including degeneracies; Fourier trace moments are elementary
%   symmetric functions of eig(V).^k, computed in O(D*P) per k.
%   Explicit determinant-minor construction is limited to N<=256;
%   no faster matrix-free action is claimed.
%   This is NOT the full D^P tensor power or a general Bogoliubov map.
%   [116] https://doi.org/10.26421/QIC5.3-3
properties={'square','unitary','normal','scalable','parameter-dependent', ...
    'eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(p),p=min(2,d);end
validateattributes(p,{'numeric'},{'scalar','integer','nonnegative','finite','<=',d});
if nargin<3||isempty(V),V=discrete_fourier(d);end
if ~isnumeric(V)||~isequal(size(V),[d d])||any(~isfinite(V(:)))
    error('fermionic_exterior_power:Base','V must be finite D-by-D.');
end
if norm(V'*V-eye(d),'fro')>2e-10*max(1,d)
    error('fermionic_exterior_power:Unitary','V must be unitary.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('fermionic_exterior_power:Mode','Invalid output mode.');
end
N=nchoosek(d,p);
if ~isfinite(N)||N>flintmax
    error('fermionic_exterior_power:Order','Dimension exceeds exact index range.');
end
switch lower(char(output))
    case 'spectral'
        z=eig(full(V));z=z./abs(z);
        U=unitary_spectral_reference('fermionic_exterior_power',d);
        U.dimension=N;
        U.finite.kind='exact';
        U.finite.eigenvalues=@()um_exterior_spectral(z,p,'eigenvalues');
        U.finite.moments=@(k)um_exterior_spectral(z,p,'moments',k);
        U.reference='Bravyi (2005), Quantum Information and Computation';
        U.notes=['Exact full P-fermion trace spectral phases from the ' ...
            'one-particle unitary; occupied modes are distinct.'];
    case {'matrix','explicit'}
        if N>256
            error('fermionic_exterior_power:Order', ...
                'Explicit minor matrix limited to order 256; use spectral.');
        end
        if p==0,U=1;return;end
        if p==d,U=det(V);return;end
        idx=nchoosek(1:d,p);
        U=zeros(N,N);
        if ~isreal(V),U=complex(U);end
        for a=1:N
            for b=1:N
                U(a,b)=det(V(idx(a,:),idx(b,:)));
            end
        end
    otherwise
        error('fermionic_exterior_power:Mode', ...
            'Only matrix and spectral modes are supported.');
end
end
