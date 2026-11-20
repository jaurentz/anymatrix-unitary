function [U,properties]=connection_cycle(N,angles,output)
%CONNECTION_CYCLE Sparse cyclic discrete-connection transport model.
%   U=CONNECTION_CYCLE(N,ANGLES) is a 2N-by-2N block permutation.
%   Its only nonzero block in block column j is R(ANGLES(j)),
%   in block row j+1 modulo N, with R=[cos -sin;sin cos].
%   ANGLES are frame-to-frame connection angles in radians; this
%   generator does NOT infer angles from an actual surface mesh.
%   Holonomy H=R(sum(ANGLES)), so its 2N eigenvalues are the Nth
%   roots of exp(+/-1i*sum(ANGLES)), known exactly.
%   Default U is sparse. 'operator' uses O(N) real rotations.
%   Geometric source for connection-angle data [94]:
%   Crane, Desbrun and Schroeder, Computer Graphics Forum 29 (2010).
%   https://doi.org/10.1111/j.1467-8659.2010.01761.x
properties={'square','unitary','normal','orthogonal','real','sparse', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(N,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(angles)
    angles=.2+.1*cos(2*pi*(1:N)'/(N+1));
end
if ~isnumeric(angles)||~isreal(angles)||~isvector(angles)|| ...
        numel(angles)~=N||any(~isfinite(angles(:)))
    error('connection_cycle:Angles','ANGLES must be N finite real angles.');
end
angles=angles(:);
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('connection_cycle:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    phi=sum(angles);k=(0:N-1)';
    z=[exp(1i*(2*pi*k+phi)/N);exp(1i*(2*pi*k-phi)/N)];
    U=unitary_spectral_reference('connection_cycle',N);
    U.dimension=2*N;U.finite.kind='exact';U.finite.eigenvalues=@()z;
    U.finite.moments=@(k)reshape(arrayfun(@(j) ...
        double(mod(j,N)==0)*cos((j/N)*phi),k(:)),size(k));
    U.reference='Discrete connection rotations; doi:10.1111/j.1467-8659.2010.01761.x';
    U.notes='Proposed cyclic block-permutation benchmark, not a published mesh operator.';
elseif any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_connection_cycle_apply(X,flag,angles);
elseif any(strcmp(output,{'matrix','explicit'}))
    j=(1:N)';d=mod(j,N)+1;c=cos(angles);s=sin(angles);
    rows=[2*d-1;2*d-1;2*d;2*d];
    cols=[2*j-1;2*j;2*j-1;2*j];
    U=sparse(rows,cols,[c;-s;s;c],2*N,2*N);
else
    error('connection_cycle:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
