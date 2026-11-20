function [U,properties]=etna_hessenberg(n,kind,seed,coupling,output)
%ETNA_HESSENBERG Structured unitary Hessenberg QR benchmark families.
%   U=ETNA_HESSENBERG(N,KIND,SEED,COUPLING) constructs paper-inspired
%   random core products for Sec. 7 of Aurentz, Mach, Vandebril and Watkins,
%   Fast and stable unitary QR algorithm, ETNA 44 (2015), 327-341.
%   https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341
%
%   KIND='I': Eq. (7.2), independently sampled modulus in (0,1) and
%     phases on the unit circle, each core first column
%     [(-1)^(k-1)*m_k*p_k; sqrt(1-m_k^2)].
%   KIND='II': real cores with first columns [cos(t_k);sin(t_k)],
%     t_k uniformly drawn from (0,pi).
%   KIND='III': N=4*B (B>=2), four copies of the *same* Type-I
%     (B-1)-core block, joined by three small-sine boundary cores.
%     A continuous unitary analog of Sec. 7 Example 7.4. Here the
%     coupling phase is fixed, rather than randomly drawn as in the
%     paper, so that the limit COUPLING->0 is exactly four copies.
%     This does not recreate the specific historical 16-eigenvalue
%     cluster at -1; it is a controlled block near-degeneracy test.
%
%   SEED defaults to 1 and is used with a temporary local RNG state;
%   the caller's RNG state is preserved. COUPLING defaults to 0.001
%   and must lie strictly between zero and one for Type III.
%   KIND may be I, II, III or 1, 2, 3. Outputs are explicit by
%   default; 'operator' gives O(N) core-factor action. Generic exact
%   eigenphases are not available.
%
%   [8] https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341
properties={'square','unitary','normal','hessenberg','scalable', ...
    'random','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(kind),kind='I';end
if isnumeric(kind) && isscalar(kind) && isfinite(kind) ...
        && kind==floor(kind) && kind>=1 && kind<=3
    names={'I','II','III'};
    kind=names{kind};
end
if ~(ischar(kind)||(isstring(kind)&&isscalar(kind)))
    error('etna_hessenberg:Kind','KIND must be I, II, or III.');
end
kind=upper(char(kind));
if ~any(strcmp(kind,{'I','II','III'}))
    error('etna_hessenberg:Kind','KIND must be I, II, or III.');
end
if strcmp(kind,'III') && (mod(n,4)~=0 || n<8)
    error('etna_hessenberg:Order','Type III requires N=4*B with B>=2.');
end
if nargin<3||isempty(seed),seed=1;end
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','finite'});
if nargin<4||isempty(coupling),coupling=0.001;end
validateattributes(coupling,{'numeric'},{'scalar','real','finite'});
if strcmp(kind,'III') && ~(coupling>0 && coupling<1)
    error('etna_hessenberg:Coupling','Type III needs 0<COUPLING<1.');
end
if strcmp(kind,'III') && 1-coupling^2==1
    error('etna_hessenberg:Coupling', ...
          'COUPLING is too small to resolve in double-precision Schur parameters.');
end
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('etna_hessenberg:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('etna_hessenberg',n);
    U.reference='Aurentz et al., ETNA 44 (2015), 327-341, Section 7';
    U.notes=['Seeded Hessenberg benchmark construction. Type III has four ' ...
        'near-decoupled identical core blocks; no generic exact trace law.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('etna_hessenberg:Mode','Supported modes: matrix, operator, spectral.');
end
prior=rng;
restore=onCleanup(@()rng(prior)); %#ok<NASGU>
rng(seed,'twister');
switch kind
    case 'I'
        alpha=random_type_i(n-1);
        beta=1;
    case 'II'
        theta=pi*rand(n-1,1);
        alpha=cos(theta); % rho=sin(theta)>0
        beta=1;
    case 'III'
        B=n/4;
        base=random_type_i(B-1);
        alpha=zeros(n-1,1);
        for j=1:4
            start=(j-1)*B;
            alpha(start+(1:B-1))=base;
        end
        alpha(B:B:n-1)=-sqrt(1-coupling^2);
        beta=-1;
end
U=opuc_matrix(n,alpha,beta,'hessenberg',[],output);
end

function alpha=random_type_i(n)
m=rand(n,1);
p=exp(2i*pi*rand(n,1));
alpha=(-1).^(0:n-1)'.*m.*conj(p);
end
