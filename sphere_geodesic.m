function [U,properties]=sphere_geodesic(n,p,q,output)
%SPHERE_GEODESIC Ambient rotation realizing sphere geodesic transport.
%   R=SPHERE_GEODESIC(N,P,Q) maps P to Q by the shortest great-circle arc
%   on S^(N-1) and restricts to Levi-Civita parallel transport between
%   tangent spaces. P,Q are real unit N-vectors, not antipodal.
%   With A=Q*P'-P*Q' and c=P'*Q, R=I+A+A^2/(1+c).
%   Spectrum: exp(+/-1i*acos(c)) and N-2 copies of 1.
%   Final OUTPUT='operator' gives O(N) forward/adjoint action;
%   'spectral' gives exact finite eigenvalues and Fourier moments.
%   Source [96], https://doi.org/10.3390/sym13112092
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(p),p=zeros(n,1);p(1)=1;end
if nargin<3||isempty(q),q=zeros(n,1);q(1)=cos(.7);q(2)=sin(.7);end
if ~isnumeric(p)||~isreal(p)||~isvector(p)||numel(p)~=n|| ...
        any(~isfinite(p(:)))||abs(norm(p)-1)>5e-12
    error('sphere_geodesic:Point','P must be a real unit N-vector.');
end
if ~isnumeric(q)||~isreal(q)||~isvector(q)||numel(q)~=n|| ...
        any(~isfinite(q(:)))||abs(norm(q)-1)>5e-12
    error('sphere_geodesic:Point','Q must be a real unit N-vector.');
end
p=p(:);q=q(:);c=max(-1,min(1,p'*q));
if 1+c<1e-12
    error('sphere_geodesic:Antipodal','The shortest geodesic is not unique near antipodes.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('sphere_geodesic:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    theta=acos(c);
    z=[exp(1i*theta);exp(-1i*theta);ones(n-2,1)];
    U=unitary_spectral_reference('sphere_geodesic',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()z;
    U.finite.moments=@(k)(n-2+2*cos(k*theta))/n;
    U.reference='Sphere Levi-Civita geodesic transport; doi:10.3390/sym13112092';
    U.notes='Ambient rotation of rank(R-I)<=2; its restriction transports tangent vectors.';
elseif any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_sphere_apply(X,flag,p,q);
elseif any(strcmp(output,{'matrix','explicit'}))
    A=q*p'-p*q';
    U=eye(n)+A+(A*A)/(1+c);
else
    error('sphere_geodesic:Mode','OUTPUT must be matrix, operator, or spectral.');
end
end
