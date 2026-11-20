function [U,properties]=dual_hahn_transform(n,alpha,beta,output)
%DUAL_HAHN_TRANSFORM Dual Hahn orthogonal polynomial transform.
%   U=DUAL_HAHN_TRANSFORM(N,ALPHA,BETA) returns an N-by-N real
%   orthogonal matrix with polynomial degree x along rows and
%   *quadratic* spectral grid lambda(k)=k*(k+ALPHA+BETA+1)
%   along columns, k=0:N-1. The dual Hahn polynomial has the
%   terminating hypergeometric representation
%   R_x(lambda(k);ALPHA,BETA,N-1)=
%     3F2(-x,-k,k+ALPHA+BETA+1;ALPHA+1,-(N-1);1).
%   By exact Hahn/dual-Hahn duality, this normalized transform
%   equals HAHN_TRANSFORM(N,ALPHA,BETA)', with signs set by the
%   positive first-row Hahn convention. Orthogonality is exact
%   in real arithmetic, but its eigenphases are not lambda(k).
%   Unlike an independent approximation using polynomial products,
%   this construction preserves the exact transpose duality.
%   Matrix formation costs O(N^3) via a symmetric Jacobi eigensolver.
%   ALPHA,BETA>-1; defaults match HAHN_TRANSFORM.
%
%   R. Koekoek, P. A. Lesky, R. F. Swarttouw (2010).
%   https://doi.org/10.1007/978-3-642-05014-5
%   Hahn/dual-Hahn duality: https://dlmf.nist.gov/18.21
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=.4;end
if nargin<3||isempty(beta),beta=.9;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(beta,{'numeric'},{'scalar','real','finite','>',-1});
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('dual_hahn_transform:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('dual_hahn_transform',n);
    U.reference='Koekoek, Lesky and Swarttouw (2010), dual Hahn class';
    U.notes=['Normalized dual Hahn transform is the adjoint of ' ...
       'the matching Hahn transform; its grid is quadratic in k. ' ...
       'The polynomial grid is not the matrix eigenvalue spectrum.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('dual_hahn_transform:Mode','Supported modes are matrix and spectral.');
end
V=um_hahn_basis(n,alpha,beta);
U=V;
end
