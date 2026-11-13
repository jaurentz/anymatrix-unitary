function [U,properties] = orthogonal_blocks(n,theta,epsilon,output)
%ORTHOGONAL_BLOCKS Real orthogonal matrix with prescribed rotation angles.
%   U=ORTHOGONAL_BLOCKS(N,THETA,EPSILON) is block diagonal with
%   floor(N/2) blocks R(t)=[cos(t),-sin(t);sin(t),cos(t)]. For odd N,
%   the last diagonal entry is EPSILON=+1 or -1 (default +1).
%   THETA must contain floor(N/2) angles (default 0.3*j).
%   Spectrum: exp(+-1i*THETA(j)), plus EPSILON if N is odd.
%   Uses: real normal EVD, conjugate pairing, spectral clusters.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [3] https://doi.org/10.13001/ela.2026.9957
%   [4] https://doi.org/10.1137/1.9781611975345
properties={'square','unitary','normal','orthogonal','real','sparse', ...
    'banded','scalable','eigenvalues','parameter-dependent'};
% Optional spectral-reference mode; preserve legacy output properties.
if nargin>=4 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('orthogonal_blocks',n,theta,epsilon);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
m=floor(n/2);
if nargin<2 || isempty(theta),theta=0.3*(1:m)';end
validateattributes(theta,{'numeric'},{'vector','real','finite'});
if numel(theta)~=m,error('orthogonal_blocks:Size','THETA must have floor(N/2) entries.');end
if nargin<3 || isempty(epsilon),epsilon=1;end
if ~(isnumeric(epsilon)&&isscalar(epsilon)&&isreal(epsilon)&&any(epsilon==[-1 1]))
    error('orthogonal_blocks:Sign','EPSILON must be +1 or -1.');
end
U=speye(n);
for k=1:m
    t=theta(k); indices=2*k-1:2*k;
    U(indices,indices)=[cos(t),-sin(t);sin(t),cos(t)];
end
if mod(n,2),U(n,n)=epsilon;end
end
