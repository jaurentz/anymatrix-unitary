function [U,properties] = orthogonal_butterfly(n,angles,output)
%ORTHOGONAL_BUTTERFLY Recursively factored real orthogonal butterfly.
%   U=ORTHOGONAL_BUTTERFLY(N,ANGLES) for power-of-two N constructs
%   the recursive Parker/Trogdon nonsimple butterfly
%      B_N = [c*I s*I;-s*I c*I] * blkdiag(B_left,B_right).
%   ANGLES contains N-1 angles in preorder: root, entire left subtree,
%   entire right subtree. Each internal node takes its own angle.
%   The default angle vector is deterministic and reproducible;
%   supply 2*pi*rand(N-1,1) for random ensembles.
%   Output='operator' gives O(N log N) forward and adjoint matvecs.
%   These are in SO(N); their law is not Haar on SO(N) in general.
%
%   T. Trogdon, Appl. Math. Lett. 2019, arXiv:1710.00087.
%   https://doi.org/10.1016/j.aml.2019.03.024
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if 2^round(log2(n))~=n
    error('orthogonal_butterfly:Order','N must be a power of two.');
end
if nargin<2||isempty(angles),angles=pi*(1:n-1)'/(n+1);end
if ~isnumeric(angles)||~isvector(angles)||numel(angles)~=n-1 ...
        ||~isreal(angles)||any(~isfinite(angles(:)))
    error('orthogonal_butterfly:Angles','ANGLES must be N-1 finite real entries.');
end
angles=angles(:);
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('orthogonal_butterfly:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('orthogonal_butterfly',n);
    U.reference='Trogdon, On spectral and numerical properties of random butterfly matrices (2019)';
    U.notes=['Deterministic angles or user-provided random angles. ' ...
        'Real orthogonal conjugate-paired eigenphases, ' ...
        'but no universal sample density.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('orthogonal_butterfly:Mode','Mode must be matrix, operator or spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_butterfly_apply(X,flag,n,angles);
else
    U=um_butterfly_apply(eye(n),'notransp',n,angles);
end
end
