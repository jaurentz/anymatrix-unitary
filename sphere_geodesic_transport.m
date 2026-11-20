function [U,properties]=sphere_geodesic_transport(n,x,y,output)
%SPHERE_GEODESIC_TRANSPORT Ambient SO(N) lift of parallel transport on S^(N-1).
%   U=SPHERE_GEODESIC_TRANSPORT(N,X,Y) uses unit vectors X and Y in
%   R^N (default e1 and cos(.8)e1+sin(.8)e2). Antipodal endpoints
%   are excluded: there is no unique shortest great-circle geodesic.
%
%   With c=X'*Y, K=Y*X'-X*Y', define the minimal planar rotation
%       U=I+K+K^2/(1+c),  c>-1.
%   This SO(N) rotation maps X to Y, fixes the perpendicular of
%   span(X,Y), and its restriction maps T_X S^(N-1) to T_Y S^(N-1)
%   by Levi-Civita parallel transport along the shortest geodesic.
%   Strictly speaking, the tangent-space map is between different
%   vector spaces; the SQUARE ambient orthogonal extension is what
%   is stored in this database.
%   The exact finite eigenvalues are e^(+/-i*acos(c)) and 1 (N-2 times);
%   the finite normalized trace moments are
%       (N-2+2*cos(k*acos(c)))/N.
%   OUTPUT='operator' applies a rank-two update in O(N) per RHS.
%   Near-antipodal pairs are rejected, since 1+c is ill-conditioned.
%
%   P.-A. Absil, R. Mahony, R. Sepulchre,
%   Optimization Algorithms on Matrix Manifolds (2008).
%   https://doi.org/10.1515/9781400830244
%   S. Fiori, Symmetry 13 (2021), 2092 (same spherical transport geometry).
%   https://doi.org/10.3390/sym13112092
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(x)
    x=zeros(n,1);x(1)=1;
end
if nargin<3||isempty(y)
    y=zeros(n,1);y(1)=cos(.8);y(2)=sin(.8);
end
if ~isnumeric(x)||~isreal(x)||~isvector(x)||numel(x)~=n ...
    ||any(~isfinite(x(:)))||abs(norm(x)-1)>1e-10
    error('sphere_geodesic_transport:Start', ...
          'X must be a real unit N-vector.');
end
if ~isnumeric(y)||~isreal(y)||~isvector(y)||numel(y)~=n ...
    ||any(~isfinite(y(:)))||abs(norm(y)-1)>1e-10
    error('sphere_geodesic_transport:End', ...
          'Y must be a real unit N-vector.');
end
x=x(:)/norm(x);y=y(:)/norm(y);
c=max(-1,min(1,x'*y));
if 1+c<=1e-12
    error('sphere_geodesic_transport:Antipodal', ...
        'Antipodal/near-antipodal points have ill-conditioned geodesic transport.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('sphere_geodesic_transport:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('sphere_geodesic_transport',n);
    theta=acos(c);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()[exp(1i*theta);exp(-1i*theta);ones(n-2,1)];
    U.finite.moments=@(k)(n-2+2*cos(k*theta))/n;
    if theta==0
        U.finite.atoms=struct('angles',0,'weights',1);
    elseif n==2
        U.finite.atoms=struct('angles',[theta;-theta], ...
            'weights',[.5;.5]);
    else
        U.finite.atoms=struct('angles',[theta;-theta;0], ...
            'weights',[1;1;n-2]/n);
    end
    U.reference='Absil, Mahony and Sepulchre, Optimization Algorithms on Matrix Manifolds (2008)';
    U.notes=['Exact finite ambient SO(n) spectrum. Tangent-space Levi-Civita ' ...
        'transport between endpoints is not an endomorphism unless the ' ...
        'tangent spaces are identified by this particular ambient lift.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('sphere_geodesic_transport:Mode', ...
          'Supported modes: matrix, operator and spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_geodesic_apply(X,flag,x,y,c);
else
    K=y*x'-x*y';
    U=eye(n)+K+(K*K)/(1+c);
end
end
