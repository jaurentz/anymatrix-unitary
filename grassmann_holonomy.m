function [H,properties]=grassmann_holonomy(n,r,frames,output)
%GRASSMANN_HOLONOMY Non-Abelian geodesic-loop tautological holonomy.
%   H=GRASSMANN_HOLONOMY(N,R,FRAMES) uses an N-by-R-by-M array
%   of orthonormal complex frames with M>=3. Each consecutive pair
%   (including the closing pair) must have nonsingular R-by-R overlap.
%   On the unique minimizing Grassmann geodesic segment from frame Y_j
%   to Y_(j+1), the canonical tautological-bundle parallel transport
%   in frame coordinates is polar(Y_(j+1)'*Y_j).
%   H is the product of these link unitaries in traversal order.
%   Under per-vertex frame gauge Y_j -> Y_j*Q_j, the closed product
%   transforms H -> Q_1'*H*Q_1, preserving eigenvalues.
%   This is NOT Levi-Civita transport on the Grassmann tangent bundle.
%   'spectral' gives an exact rank-one Bargmann phase, but labels
%   generic non-Abelian rank-R>=2 spectra unknown analytically.
%   [97] Mukunda and Simon, Ann. Phys. 228 (1993), 205-268.
%   https://doi.org/10.1006/aphy.1993.1093
%   [98] Bendokat, Zimmermann, Absil, Adv. Comput. Math. 50 (2024).
%   https://doi.org/10.1007/s10444-023-10090-8
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,H=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2||isempty(r),r=min(2,n-1);end
validateattributes(r,{'numeric'},{'scalar','integer','positive','finite'});
if r>=n,error('grassmann_holonomy:Order','R must be smaller than N.');end
if nargin<3||isempty(frames)
    frames=zeros(n,r,3);
    frames(:,:,1)=eye(n,r);
    frames(:,:,2)=frames(:,:,1);
    a=.55;frames(1,1,2)=cos(a);frames(r+1,1,2)=sin(a);
    frames(:,:,3)=frames(:,:,1);
    b=.7;
    if r==1
        frames(1,1,3)=cos(b);frames(2,1,3)=exp(.6i)*sin(b);
    else
        frames(2,2,3)=cos(b);frames(r+1,2,3)=exp(.6i)*sin(b);
    end
end
if ~isnumeric(frames)||~isfloat(frames)|| ...
        size(frames,1)~=n||size(frames,2)~=r||size(frames,3)<3 ...
        ||any(~isfinite(frames(:)))
    error('grassmann_holonomy:Frames','FRAMES must be finite N-by-R-by-M, M>=3.');
end
m=size(frames,3);
for j=1:m
    Q=frames(:,:,j);
    if norm(Q'*Q-eye(r),'fro')>1e-10*max(1,r)
        error('grassmann_holonomy:Frames','Every frame must have orthonormal columns.');
    end
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('grassmann_holonomy:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','explicit','spectral'}))
    error('grassmann_holonomy:Mode','OUTPUT must be matrix or spectral.');
end
P=eye(r);
for j=1:m
    Y=frames(:,:,j);
    Z=frames(:,:,mod(j,m)+1);
    [L,S,V]=svd(Z'*Y,'econ');
    if min(diag(S))<=1e-12
        error('grassmann_holonomy:CutLocus', ...
            'Adjacent subspaces have a singular or nearly singular overlap.');
    end
    P=(L*V')*P;
end
if strcmp(output,'spectral')
    H=unitary_spectral_reference('grassmann_holonomy',r);
    H.reference='Mukunda-Simon (1993), doi:10.1006/aphy.1993.1093; Bendokat et al. (2024), doi:10.1007/s10444-023-10090-8';
    H.notes='Tautological U(R) holonomy of piecewise-minimizing Grassmann geodesics; spectra are gauge invariant.';
    if r==1
        H.finite.kind='exact';z=P;
        H.finite.eigenvalues=@()z;
        H.finite.moments=@(k)z.^k;
        H.notes='Exact rank-one Bargmann/Pancharatnam phase from normalized overlaps.';
    else
        H.notes=[H.notes ' No generic closed-form non-Abelian spectrum is supplied.'];
    end
else
    H=P;
end
end
