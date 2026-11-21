function [U,properties]=grassmann_tangent_holonomy(n,r,frames,field,mode)
%GRASSMANN_TANGENT_HOLONOMY Intrinsic Grassmann Levi-Civita holonomy.
%   U=GRASSMANN_TANGENT_HOLONOMY(N,R,FRAMES,FIELD,MODE)
%   constructs the intrinsic tangent parallel transport around a
%   geodesic polygon in real Gr(R,N) or complex Gr(R,N), not the
%   N-by-N ambient lift GRASSMANN_GEODESIC nor the tautological
%   R-by-R Berry/Wilczek-Zee holonomy GRASSMANN_HOLONOMY.
%
%   FRAMES is N-by-R-by-M, M>=3, with orthonormal real/complex
%   columns Y_j; closing geodesic from last to first is implicit.
%   Every consecutive overlap Y_(j+1)'*Y_j must be invertible so
%   the principal-angle geodesic lies strictly within the cut locus.
%   Default N=3,R=1 and complex three-vertex loop.
%   FIELD='real' or 'complex'; for omitted FIELD, real user-supplied
%   frames select real and otherwise complex. N=2..16, 1<=R<N.
%
%   Horizontal tangent at Y is V=Qperp*B, Y'*V=0, with Qperp=null(Y')
%   and B of order (N-R)-by-R. Its Frobenius metric is
%       <V,W>=trace(B'*C).
%   For each minimizing geodesic Y->Z, align Z to Y via
%       Ralign=polar(Z'*Y), Zend=Z*Ralign,
%   form the horizontal log Delta and skew(-Hermitian) generator
%       Omega=Delta*Y'-Y*Delta',  A=exp(Omega), A*Y=Zend.
%   Along Y(t)=exp(t*Omega)*Y, the Levi-Civita parallel transport
%   of a horizontal tangent is V(t)=exp(t*Omega)*V.
%   In FIXED orthonormal coordinate frames the intrinsically
%   unitary tangent link is
%       B_next=S*B*Ralign', S=Qperp_end'*A*Qperp_start,
%       vec(B_next)=kron(conj(Ralign),S)*vec(B).
%   The returned (N-R)*R square complex unitary (or real orthogonal)
%   is the PRODUCT of those links in geodesic polygon order.
%   The output matrix order is tangent dimension, NOT N and NOT R.
%
%   MODE='matrix' (default): intrinsic square tangent holonomy.
%   MODE='model': per-leg ambient unitary A, horizontal velocity,
%                 gauge/alignment R, complement S, tangent links,
%                 and total tangent holonomy.
%   MODE='spectral': generic real/complex Grassmann analytic finite
%                    trace eigenphase law UNKNOWN; complex CP^1
%                    has exact scalar tangent holonomy
%                    conj(Berry_phase)^2, not Berry_phase itself.
%
%   In CP^1=Gr_C(1,2), Fubini-Study Gauss curvature is 4;
%   its tangent holonomy is the SQUARE-INVERSE of tautological
%   line Berry transport (in these sign conventions). For CP^m,
%   m>=2, the tangent U(m) holonomy can be non-Abelian.
%
%   [97] Edelman, Arias, Smith, SIAM JMAA 20 (1998), 303-353.
%        https://doi.org/10.1137/S0895479895290954
%   [111] Bendokat, Zimmermann, Absil (2024), Grassmann handbook.
%         https://doi.org/10.1007/s10444-023-10090-8
%   [110] Mukunda and Simon (1993), geometrical phase for contrast.
%         https://doi.org/10.1006/aphy.1993.1093
properties={'square','unitary','normal','parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
       &&n==floor(n)&&n>=2&&n<=16)
    error('grassmann_tangent_holonomy:Order','N must be 2..16.');
end
if nargin<2||isempty(r),r=min(2,n-1);end
if ~(isnumeric(r)&&isscalar(r)&&isreal(r)&&isfinite(r) ...
        &&r==floor(r)&&r>=1&&r<n)
    error('grassmann_tangent_holonomy:Rank','R must be in 1..N-1.');
end
if nargin<3||isempty(frames)
    frames=zeros(n,r,3);
    frames(:,:,1)=eye(n,r);
    frames(:,:,2)=eye(n,r);
    frames(1,1,2)=cos(.43);
    frames(r+1,1,2)=sin(.43);
    frames(:,:,3)=eye(n,r);
    if r==1
        frames(1,1,3)=cos(.51);
        frames(r+1,1,3)=exp(.7i)*sin(.51);
    else
        frames(2,2,3)=cos(.51);
        frames(r+1,2,3)=exp(.7i)*sin(.51);
    end
end
if nargin<4||isempty(field)
    if isreal(frames),field='real';else,field='complex';end
end
if nargin<5||isempty(mode),mode='matrix';end
if ~(ischar(field)||(isstring(field)&&isscalar(field))) ...
    || ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('grassmann_tangent_holonomy:Mode','FIELD and MODE must be text.');
end
field=lower(char(field));mode=lower(char(mode));
if ~any(strcmp(field,{'real','complex'}))
    error('grassmann_tangent_holonomy:Field','FIELD is real or complex.');
end
if ~any(strcmp(mode,{'matrix','model','spectral'}))
    error('grassmann_tangent_holonomy:Mode','Unknown output mode.');
end
if ~(isnumeric(frames)&&isfloat(frames)&&ndims(frames)<=3 ...
    &&size(frames,1)==n&&size(frames,2)==r ...
    &&size(frames,3)>=3&&size(frames,3)<=20 ...
    &&all(isfinite(frames(:))))
    error('grassmann_tangent_holonomy:Frames', ...
        'FRAMES must be finite N-by-R-by-M array, 3<=M<=20.');
end
if strcmp(field,'real')&&~isreal(frames)
    error('grassmann_tangent_holonomy:Frames', ...
        'Real Grassmann tangent transport requires real frames.');
end
m=size(frames,3);dim=r*(n-r);
for i=1:m
    Y=frames(:,:,i);Z=frames(:,:,mod(i,m)+1);
    if norm(Y'*Y-eye(r),'fro')>1e-10*max(1,r) ...
        || min(svd(Z'*Y))<=1e-10
        error('grassmann_tangent_holonomy:CutLocus', ...
            'Orthonormal frames and nonsingular adjacent overlap required.');
    end
end
if strcmp(mode,'spectral')
    blank=struct('kind','unknown','density',[],'support',[], ...
        'atoms',[],'moments',[]);
    finite=struct('kind','unknown','eigenvalues',[], ...
        'moments',[],'atoms',[]);
    notes=['Analytic finite tangent holonomy eigenphase laws are ', ...
      'not generally known for Grassmann geodesic polygons. ', ...
      'Tangent holonomy differs from the ambient n-by-n geodesic ', ...
      'lift and from rank-r tautological Berry/Wilczek-Zee holonomy.'];
    if n==2 && r==1 && strcmp(field,'complex')
        berry=1;
        for i=1:m
            overlap=frames(:,mod(i,m)+1)'*frames(:,i);
            berry=berry*overlap/abs(overlap);
        end
        lambda=conj(berry)^2;
        finite.kind='exact';
        finite.eigenvalues=@()lambda;
        finite.moments=@(k)lambda.^k;
        finite.atoms=struct('angles',angle(lambda),'weights',1);
        notes=['Exact complex CP^1 tangent scalar: inverse-square ', ...
            'of tautological Berry line phase, reflecting ', ...
            'Fubini-Study curvature 4 (rather than line curvature 2).'];
    end
    U=struct('family','grassmann_tangent_holonomy', ...
        'dimension',dim, ...
        'normalization','probability (density with respect to dtheta)', ...
        'finite',finite,'limit',blank,'ensemble',blank, ...
        'underlying',blank, ...
        'reference','Edelman, Arias, Smith (1998); Bendokat, Zimmermann, Absil (2024)', ...
        'notes',notes,'ambient_dimension',n, ...
        'rank',r,'field',field);
    return
end
U=eye(dim);
if strcmp(mode,'model')
    links=zeros(dim,dim,m);
    lifts=zeros(n,n,m);
    alignments=zeros(r,r,m);
    complements=zeros(n-r,n-r,m);
    velocities=zeros(n,r,m);
end
for i=1:m
    Y=frames(:,:,i);Z=frames(:,:,mod(i,m)+1);
    [T,A,Ralign,S,Delta]=um_grassmann_tangent_leg(Y,Z);
    U=T*U;
    if strcmp(mode,'model')
        links(:,:,i)=T;
        lifts(:,:,i)=A;
        alignments(:,:,i)=Ralign;
        complements(:,:,i)=S;
        velocities(:,:,i)=Delta;
    end
end
if strcmp(field,'real'),U=real(U);end
if strcmp(mode,'model')
    U=struct('U',U,'links',links,'ambient_lifts',lifts, ...
        'alignments',alignments,'complement_maps',complements, ...
        'geodesic_initial_velocities',velocities, ...
        'frames',frames,'rank',r,'field',field, ...
        'tangent_dimension',dim, ...
        'unitarity_residual',norm(U'*U-eye(dim),'fro'), ...
        'notes',['Intrinsic Levi-Civita tangent holonomy acting ', ...
          'on vec(B) where B is (N-R)-by-R. Ambient lifts ', ...
          'and rank-R tautological transports are different matrices.']);
end
end
