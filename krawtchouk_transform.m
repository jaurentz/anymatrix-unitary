function [U,properties]=krawtchouk_transform(n,p,output)
%KRAWTCHOUK_TRANSFORM Real orthogonal binomial Krawtchouk transform.
%   U=KRAWTCHOUK_TRANSFORM(N,P) is the N-by-N normalized Krawtchouk
%   transform for a binomial measure with N-1 trials and success
%   probability 0<P<1 (default P=1/2). It is a real symmetric
%   orthogonal involution U*U=I.
%
%   Construct a Jacobi matrix with diagonal
%     J(x,x)=(N-1-x)*P+x*(1-P), x=0,...,N-1,
%   and negative off-diagonal
%     J(x,x+1)=-sqrt(P*(1-P)*(x+1)*(N-1-x)).
%   Its eigenvalues are exactly 0,...,N-1. Choose each eigenvector
%   to have positive FIRST COMPONENT; the resulting matrix U contains
%   the normalized Krawtchouk functions, with the polynomial-degree
%   index along rows and the grid index along columns in the
%   self-dual convention. These eigenvalues are the grid nodes of J,
%   NOT the eigenvalues of the Krawtchouk transform U (which are +/-1).
%
%   The eigensolver-based construction costs O(N^3) and does NOT claim
%   a fast Krawtchouk algorithm. Very high degrees or probabilities
%   close to 0 or 1 can make sign-setting ill-conditioned; reject
%   cases where the computed first-row sign becomes unresolved.
%
%   P. Feinsilver and R. Schott, Krawtchouk transforms and convolutions.
%   https://doi.org/10.1007/s13373-018-0132-2
properties={'square','unitary','normal','orthogonal','real','symmetric', ...
            'involutory','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(p),p=.5;end
validateattributes(p,{'numeric'},{'scalar','real','finite','>',0,'<',1});
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('krawtchouk_transform:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('krawtchouk_transform',n);
    U.reference='Feinsilver and Schott, Krawtchouk transforms and convolutions (2018)';
    U.finite.kind='exact';
    U.finite.eigenvalues=@()[ones(ceil(n/2),1);-ones(floor(n/2),1)];
    U.finite.atoms=struct('angles',[0;pi], ...
        'weights',[ceil(n/2);floor(n/2)]/n);
    U.finite.moments=@(k)reshape(arrayfun(@(j) ...
        (ceil(n/2)+(-1).^j*floor(n/2))/n,k(:)),size(k));
    U.notes=['U is a symmetric orthogonal involution, with +/-1 finite ' ...
        'eigenvalues. The Jacobi nodes 0:N-1 are NOT eigenphases of U.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('krawtchouk_transform:Mode','Only matrix or spectral supported.');
end
x=(0:n-1)';
diagonal=p*(n-1-x)+(1-p)*x;
off=-sqrt(p*(1-p)*(1:n-1)'.*(n-(1:n-1)'));
J=diag(diagonal)+diag(off,1)+diag(off,-1);
[V,L]=eig(J,'vector');
[~,idx]=sort(real(L));
V=real(V(:,idx));
s=sign(V(1,:));
if any(s==0)
    error('krawtchouk_transform:Precision', ...
        'First components of the Jacobi eigenbasis are unresolved.');
end
U=bsxfun(@times,V,s);
end
