function [U,properties]=hyperbolic_triangle_holonomy(n,vertices,output)
%HYPERBOLIC_TRIANGLE_HOLONOMY H^2 Levi-Civita holonomy via Gauss-Bonnet.
%   U=HYPERBOLIC_TRIANGLE_HOLONOMY(N,VERTICES) constructs the SO(2)
%   tangent-space parallel-transport holonomy of the geodesic triangle
%   with vertices Z1,Z2,Z3 in the Poincare disk |Z|<1. VERTICES
%   is a complex 3-vector, default [0;.45;.35i]. The curvature is -1.
%   With P=prod_j(1-conj(Z_j)*Z_(j+1)), Z4=Z1, oriented area
%       A = -2*angle(P),  -pi<A<pi.
%   Gauss-Bonnet then gives the holonomy rotation R(-A), where
%       R(theta)=[cos(theta) -sin(theta);sin(theta) cos(theta)].
%   OUTPUT is R(-A) padded by I_(N-2), N>=2. Padding is a
%   *benchmark direct sum*, not a claim that the hyperbolic
%   tangent space is N dimensional. In particular N=2 is geometric.
%   Eigenvalues exp(+-i*A), and N-2 additional unit eigenvalues.
%   N-by-N output and matrix-free O(N) action; the latter consists
%   of the leading 2x2 rotation and identity spectators.
%
%   T. Needham, Visual Differential Geometry and Forms (2021),
%   Gauss-Bonnet and holonomy; https://doi.org/10.1515/9780691219899
%   J. Voight, Quaternion Algebras (2021), hyperbolic triangle area,
%   Section 33.6; https://doi.org/10.1007/978-3-030-56694-4
properties={'square','unitary','normal','orthogonal','real', ...
            'scalable','parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(vertices),vertices=[0;.45;.35i];end
if ~isnumeric(vertices)||~isvector(vertices)||numel(vertices)~=3 ...
        ||any(~isfinite(vertices(:)))||any(abs(vertices(:))>=1)
    error('hyperbolic_triangle_holonomy:Vertices', ...
          'VERTICES must be three finite complex points in the open disk.');
end
vertices=vertices(:);
vnext=vertices([2 3 1]);
P=prod(1-conj(vertices).*vnext);
area=-2*angle(P);
theta=-area;
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('hyperbolic_triangle_holonomy:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('hyperbolic_triangle_holonomy',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()[exp(1i*theta);exp(-1i*theta);ones(n-2,1)];
    U.finite.moments=@(k)(n-2+2*cos(k*theta))/n;
    if abs(theta)<5e-14
        U.finite.atoms=struct('angles',0,'weights',1);
    elseif n==2
        U.finite.atoms=struct('angles',[theta;-theta], ...
                              'weights',[.5;.5]);
    else
        U.finite.atoms=struct('angles',[theta;-theta;0], ...
                              'weights',[1;1;n-2]/n);
    end
    U.reference='Needham, Visual Differential Geometry and Forms (2021)';
    U.notes=['Signed hyperbolic area equals minus the tangent ' ...
        'rotation angle (curvature -1). Only the leading SO(2) ' ...
        'block is intrinsic; remaining identity dimensions are ' ...
        'explicit benchmark padding.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('hyperbolic_triangle_holonomy:Mode', ...
          'Supported modes: matrix, operator, spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_hyperbolic_holonomy_apply(X,flag,n,theta);
else
    U=eye(n);
    U(1:2,1:2)=[cos(theta) -sin(theta);sin(theta) cos(theta)];
end
end
