function [G,properties]=grassmann_geodesic(n,angles,frames,output)
%GRASSMANN_GEODESIC Ambient unitary lift of a Grassmann geodesic.
%   G=GRASSMANN_GEODESIC(N,ANGLES,FRAMES) uses r=length(ANGLES),
%   2*r<=N, and N-by-2r orthonormal columns FRAMES=[Y,U].
%   The tangent H=U*diag(ANGLES) defines the Grassmann geodesic
%       Y(t)=Y*cos(t*diag(ANGLES))+U*sin(t*diag(ANGLES)).
%   G at t=1 is the ambient unitary acting by r planar rotations
%   on [Y,U] and identically on its complement. It is NOT the
%   parallel-transport map on the Grassmann tangent bundle.
%   Eigenvalues: exp(+/-1i*ANGLES), plus N-2r unit eigenvalues.
%   Complex FRAMES yield complex unitary G; real FRAMES yield orthogonal G.
%   'operator': O(Nr) low-rank forward/adjoint; 'spectral': exact law.
%   Edelman, Arias & Smith, SIAM JMAA 20 (1998), 303-353 [93].
%   https://doi.org/10.1137/S0895479895290954
properties={'square','unitary','normal','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,G=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(angles)
    angles=.25+.2*(1:min(3,floor(n/2)))';
end
if ~isnumeric(angles)||~isreal(angles)||~isvector(angles)|| ...
        isempty(angles)||any(~isfinite(angles(:)))||2*numel(angles)>n
    error('grassmann_geodesic:Angles','ANGLES must be r real finite angles, 1<=2r<=N.');
end
angles=angles(:);r=numel(angles);
if nargin<3||isempty(frames),frames=eye(n,2*r);end
if ~isnumeric(frames)||ndims(frames)>2|| ...
        ~isequal(size(frames),[n 2*r])||any(~isfinite(frames(:)))|| ...
        norm(frames'*frames-eye(2*r),'fro')>1e-10*max(1,r)
    error('grassmann_geodesic:Frames','FRAMES must have 2r orthonormal N-vectors.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('grassmann_geodesic:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    z=[exp(1i*angles);exp(-1i*angles);ones(n-2*r,1)];
    G=unitary_spectral_reference('grassmann_geodesic',n);
    G.finite.kind='exact';G.finite.eigenvalues=@()z;
    G.finite.moments=@(k)reshape(arrayfun(@(j) ...
        (n-2*r+2*sum(cos(j*angles)))/n,k(:)),size(k));
    G.reference='Edelman, Arias, Smith, SIAM JMAA 20 (1998); doi:10.1137/S0895479895290954';
    G.notes='Ambient unitary geodesic lift, not intrinsic tangent parallel transport.';
elseif any(strcmp(output,{'operator','implicit'}))
    G=@(X,flag)um_grassmann_apply(X,flag,frames,angles);
elseif any(strcmp(output,{'matrix','explicit'}))
    c=diag(cos(angles));s=diag(sin(angles));
    C=[c -s;s c];
    G=eye(n)+frames*(C-eye(2*r))*frames';
else
    error('grassmann_geodesic:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
