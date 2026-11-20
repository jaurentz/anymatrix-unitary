function [U,properties]=unitary_rank_one_update(n,phi,alpha,beta,v,output)
%UNITARY_RANK_ONE_UPDATE Multiplicative rank-one unitary perturbation.
%   U=UNITARY_RANK_ONE_UPDATE(N,PHI,ALPHA,BETA,V) constructs
%       H*(I+(exp(1i*PHI)-1)*v*v'), norm(v)=1,
%   where H is the N-by-N unitary upper Hessenberg matrix specified
%   by N-1 Verblunsky/core coefficients ALPHA and boundary BETA.
%   The default v=e_N preserves Hessenberg structure by modifying only
%   the final column. Arbitrary v is permitted but need not preserve it.
%   Both U and H are exactly unitary in exact arithmetic; the matrix
%   difference has rank <= 1, and det(U)/det(H)=exp(1i*PHI).
%   When ALPHA=0 and v=e_N, roots have explicit cyclic phases.
%   'operator' returns a core-factor plus rank-one action in O(N).
%
%   P. J. Forrester and E. M. Rains, Jacobians and rank 1 perturbations
%   relating to unitary Hessenberg matrices, IMRN 2006, article 48306.
%   https://doi.org/10.1155/IMRN/2006/48306
%   B. Simon, Orthogonal Polynomials on the Unit Circle, 2005.
%   https://doi.org/10.1090/coll054.1
%   https://doi.org/10.1090/coll/054.2
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(phi),phi=.7;end
validateattributes(phi,{'numeric'},{'scalar','real','finite'});
if nargin<3||isempty(alpha),alpha=.2*exp(.3i*(1:n-1)');end
if ~isnumeric(alpha)||~isvector(alpha)||numel(alpha)~=n-1 ...
        ||any(~isfinite(alpha(:)))||any(abs(alpha(:))>=1)
    error('unitary_rank_one_update:Alpha','ALPHA must be N-1 strict-disk coefficients.');
end
alpha=alpha(:);
if nargin<4||isempty(beta),beta=1;end
if ~isnumeric(beta)||~isscalar(beta)||~isfinite(beta) ...
        ||abs(abs(beta)-1)>1e-12
    error('unitary_rank_one_update:Beta','BETA must have unit modulus.');
end
if nargin<5||isempty(v),v=zeros(n,1);v(n)=1;end
if ~isnumeric(v)||~isvector(v)||numel(v)~=n ...
        ||any(~isfinite(v(:)))||norm(v(:))==0
    error('unitary_rank_one_update:Vector','V must be a nonzero finite N-vector.');
end
v=v(:)/norm(v);
if nargin<6||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('unitary_rank_one_update:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('unitary_rank_one_update',n);
    U.reference='Forrester and Rains, IMRN 2006, doi:10.1155/IMRN/2006/48306';
    if all(alpha==0) && norm(v-[zeros(n-1,1);1])<1e-14
        z=exp(1i*(phi-angle(beta)+2*pi*(0:n-1)')/n);
        U.finite.kind='exact';
        U.finite.eigenvalues=@()z;
        U.finite.moments=@(k)reshape(arrayfun(@(j)mean(z.^j),k(:)),size(k));
        U.notes='Exact cyclic finite eigenphase law for zero interior Schur coefficients.';
    else
        U.notes=['Multiplicative rank-one perturbation preserves unitarity, ' ...
            'but no general finite eigenphase counting formula is asserted.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('unitary_rank_one_update:Mode','Supported modes: matrix, operator, spectral.');
end
tau=exp(1i*phi)-1;
if any(strcmp(output,{'operator','implicit'}))
    H=opuc_matrix(n,alpha,beta,'hessenberg',[],'operator');
    U=@(X,flag)um_rank_one_apply(X,flag,H,n,v,tau);
else
    H=opuc_matrix(n,alpha,beta,'hessenberg');
    U=H+tau*(H*v)*v';
end
end
