function [U,properties]=hahn_transform(n,alpha,beta,output)
%HAHN_TRANSFORM Normalized discrete Hahn polynomial transform.
%   U=HAHN_TRANSFORM(N,ALPHA,BETA) returns the real N-by-N
%   orthogonal Hahn transform. It has polynomial degree k=0:N-1
%   along ROWS, and grid x=0:N-1 along COLUMNS.
%   ALPHA,BETA>-1 (defaults 0.4,0.9). We form a symmetric
%   Jacobi difference matrix J of order N:
%      B(x)=(x+ALPHA+1)*(x-(N-1)),
%      D(x)=x*(x-BETA-N),
%      J(x,x)=-B(x)-D(x),
%      J(x,x+1)=-sqrt(B(x)*D(x+1)).
%   Its eigenvalues are lambda_k=k*(k+ALPHA+BETA+1).
%   Column k of the orthogonal eigenvector matrix V consists
%   of weighted/normalized Hahn polynomial values on x=0:N-1.
%   Fix V(1,k)>0, as Q_k(0)=1. U=V' is the analysis transform.
%   The lambda_k are eigenvalues of J, NOT generally of U.
%   Matrix formation uses dense EIG, O(N^3); a fast algorithm
%   and an analytic phase density are not asserted.
%   DUAL_HAHN_TRANSFORM(N,ALPHA,BETA) returns U', exactly
%   reflecting the Hahn / dual-Hahn polynomial duality.
%
%   R. Koekoek, P. A. Lesky, R. F. Swarttouw,
%   Hypergeometric Orthogonal Polynomials and Their q-Analogues (2010).
%   https://doi.org/10.1007/978-3-642-05014-5
%   NIST Digital Library of Mathematical Functions, Hahn class:
%   https://dlmf.nist.gov/18.22
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
    error('hahn_transform:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('hahn_transform',n);
    U.reference='Koekoek, Lesky and Swarttouw (2010), Hahn class';
    U.notes=['A normalized real orthogonal polynomial transform; ' ...
       'Jacobi eigenvalues k(k+alpha+beta+1) are polynomial nodes, ' ...
       'NOT eigenvalues of the generally nonsymmetric transform.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('hahn_transform:Mode','Supported modes are matrix and spectral.');
end
V=um_hahn_basis(n,alpha,beta);
U=V';
end
