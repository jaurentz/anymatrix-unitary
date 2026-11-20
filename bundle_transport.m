function [U,properties]=bundle_transport(n,permutation,links,output)
%BUNDLE_TRANSPORT Sparse unitary block-permutation with fiber transports.
%   U=BUNDLE_TRANSPORT(N,PERM,LINKS) maps e_j tensor v to
%   e_PERM(j) tensor LINKS(:,:,j)*v, where PERM permutes 1:N.
%   LINKS contains N unitary d-by-d matrices. Order is N*d.
%   Default: one cycle and deterministic 2-by-2 real rotations.
%   'operator' applies forward/adjoint in O(N*d^2) per RHS.
%   'spectral' gives exact eigenphases/moments from small holonomies.
%   See https://doi.org/10.1002/cpa.21395 for connection background.
properties={'square','unitary','normal','sparse','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(permutation),permutation=[2:n 1];end
if ~(isnumeric(permutation)&&isvector(permutation) ...
        &&numel(permutation)==n&&isreal(permutation) ...
        &&all(isfinite(permutation(:))) ...
        &&isequal(sort(double(permutation(:))),(1:n)'))
    error('bundle_transport:Permutation','PERM must permute 1:N.');
end
permutation=double(permutation(:));
if nargin<3||isempty(links)
    links=zeros(2,2,n);
    for j=1:n
        a=.2+.1*cos(2*pi*j/(n+1));
        links(:,:,j)=[cos(a),-sin(a);sin(a),cos(a)];
    end
end
d=size(links,1);
if ~isnumeric(links)||isempty(links)||d<1 ...
        ||size(links,2)~=d||size(links,3)~=n ...
        ||any(~isfinite(links(:)))
    error('bundle_transport:Links','LINKS must be finite d-by-d-by-N.');
end
for j=1:n
    R=links(:,:,j);
    if norm(R'*R-eye(d),'fro')>2e-10*max(1,d)
        error('bundle_transport:Unitary','Every link must be unitary.');
    end
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('bundle_transport:Mode','OUTPUT must be a scalar string.');
end
switch lower(char(output))
    case 'spectral'
        U=unitary_spectral_reference('bundle_transport',n);
        U.dimension=n*d;
        U.finite.kind='exact';
        U.finite.eigenvalues=@()um_bundle_spectral(permutation,links,'eigenvalues');
        U.finite.moments=@(k)um_bundle_spectral(permutation,links,'moments',k);
        U.reference='Block-monomial cycle reduction; Singer and Wu (2012)';
        U.notes=['Exact trace eigenphases from supplied small fiber holonomies; ' ...
            'not a mesh-derived Levi-Civita transport.'];
    case {'operator','implicit'}
        U=@(X,flag)um_bundle_transport_apply(X,flag,permutation,links);
    case {'matrix','explicit'}
        U=spalloc(n*d,n*d,n*d*d);
        for j=1:n
            to=(permutation(j)-1)*d+(1:d);
            from=(j-1)*d+(1:d);
            U(to,from)=links(:,:,j);
        end
    otherwise
        error('bundle_transport:Mode','Use matrix, operator or spectral.');
end
end
