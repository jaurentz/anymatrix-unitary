function [U,properties]=lattice_wilson_loop(n,links,output)
%LATTICE_WILSON_LOOP Matrix-valued non-Abelian lattice-gauge holonomy.
%   W=LATTICE_WILSON_LOOP(N,LINKS) returns the ordered product
%       W=L1*L2*...*Lm
%   of oriented unitary lattice gauge links on a closed m-edge contour.
%   LINKS is N-by-N-by-M and must contain unitary matrices. Default
%   LINKS are four deterministic noncommuting SU(N) matrices generated
%   as expm(i*H_j), with traceless Hermitian H_j. N=1 is trivial.
%   An edge reversed in orientation uses L_j' instead of L_j.
%
%   The lattice gauge transformation for edge j from site j to j+1
%   is L_j -> G_j * L_j * G_(j+1)', G_(m+1)=G_1. The Wilson loop
%   transforms as W -> G_1*W*G_1': eigenvalues and tr(W) are gauge
%   invariant, whereas individual matrix entries generally are not.
%   For the supplied unitary links W is exactly unitary; its determinant
%   equals prod(det(L_j)); the default traceless-generator links lie
%   in SU(N) (up to floating-point errors).
%
%   OUTPUT='operator' exposes the ordered forward/adjoint factors
%   in O(M*N^2) work per RHS, without forming the product W.
%   OUTPUT='links' exposes the N-by-N-by-M link array for checking
%   gauge covariance. 'spectral' records unknown generic finite
%   eigenphases. The observable tr(W)/N is gauge invariant but
%   is not a normalized probability density on the unit circle.
%
%   A. Ipp, D. Müller, Implicit schemes for real-time lattice
%   gauge theory, Eur. Phys. J. C 78 (2018), 884.
%   https://doi.org/10.1140/epjc/s10052-018-6323-x
properties={'square','unitary','normal','complex','scalable', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin>=3 && (ischar(output)||(isstring(output)&&isscalar(output))) ...
        && strcmpi(char(output),'spectral') && (nargin<2||isempty(links))
    U=unitary_spectral_reference('lattice_wilson_loop',n);
    U.reference='Ipp and Muller, Eur. Phys. J. C 78 (2018), 884';
    U.notes=['Default four-link SU(n) Wilson loop. Eigenphases ' ...
        'are gauge-conjugacy invariants but have no generic exact ' ...
        'finite counting law. Spectral metadata is independent of ' ...
        'costly matrix-exponential link construction.'];
    return
end
if nargin<2||isempty(links)
    links=default_links(n);
else
    if ~isnumeric(links)||size(links,1)~=n||size(links,2)~=n ...
            ||ndims(links)>3||size(links,3)<2||any(~isfinite(links(:)))
        error('lattice_wilson_loop:Links', ...
              'LINKS must be finite N-by-N-by-M with M>=2.');
    end
    for j=1:size(links,3)
        L=links(:,:,j);
        if norm(L'*L-eye(n),'fro')>1e-10*max(1,sqrt(n))
            error('lattice_wilson_loop:Unitary', ...
                  'Every link must be unitary.');
        end
    end
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('lattice_wilson_loop:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('lattice_wilson_loop',n);
    U.reference='Ipp and Muller, Eur. Phys. J. C 78 (2018), 884';
    U.notes=['Gauge links may be noncommuting; eigenphases depend on ' ...
        'the closed contour and link matrices. Only conjugacy ' ...
        'class and traces are basepoint gauge invariant. Default ' ...
        'links lie in SU(n), user links need only be unitary.'];
    return
end
switch output
    case 'links'
        U=links;
    case {'operator','implicit'}
        U=@(X,flag)um_wilson_apply(X,flag,links);
    case {'matrix','explicit'}
        U=eye(n);
        for j=1:size(links,3)
            U=U*links(:,:,j);
        end
    otherwise
        error('lattice_wilson_loop:Mode', ...
              'Supported modes: matrix, operator, links, spectral.');
end
end

function links=default_links(n)
% Deterministic unitary plaquette links, no RNG state changes.
m=4;links=complex(zeros(n,n,m));
d=(0:n-1)'-(n-1)/2;
A=diag(d/max(1,n-1));
S=zeros(n);
for k=1:n-1,S(k,k+1)=1;end
B=(S+S')/2;C=1i*(S-S')/2;
for j=1:m
    H=(.3+.17*j)*A + (.5*cos(.7*j))*B + (.6*sin(.5*j))*C;
    H=(H+H')/2;
    H=H-(trace(H)/n)*eye(n);
    links(:,:,j)=expm(1i*H);
end
end
