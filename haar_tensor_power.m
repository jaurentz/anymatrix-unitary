function [U,properties] = haar_tensor_power(d,t,seed,output)
%HAAR_TENSOR_POWER The same Haar U(d) repeated on t tensor factors.
%   U=HAAR_TENSOR_POWER(D,T,SEED) returns a dense Kronecker power.
%   U=HAAR_TENSOR_POWER(D,T,SEED,'operator') returns a matrix-free
%   action using O(T*d^(T+1)) arithmetic per RHS and O(d^T) work space,
%   without ever materializing the order-d^T Kronecker matrix.
%   U=HAAR_TENSOR_POWER(D,T,SEED,'spectral') requires a reproducible SEED;
%   exact finite sample moments and eigenvalues are computed from a
%   small d-by-d base unitary, not from the enormous tensor power.
%   Normalized moments: (trace(V^k)/d)^t. Full tensor, not Sym^t(V).
%   Ref: L. Marcinkowska et al., arXiv:2605.31168 (2026).
%   https://arxiv.org/abs/2605.31168
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(t),t=2;end
validateattributes(t,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3,seed=[];end
if nargin<4||isempty(output),output='matrix';end
n=d^t;
if ~isfinite(n)||n>flintmax
    error('haar_tensor_power:Order','d^t exceeds exact index range.');
end
if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('haar_tensor_power',d,t,seed);
    return
end
implicit=um_is_operator(output);
if ~implicit && n>4096
    error('haar_tensor_power:Order','Explicit order >4096; use operator mode.');
end
V=haar_unitary(d,seed);
if implicit
    factors=repmat({V},1,t);
    U=@(X,flag)um_tensor_apply(X,flag,factors);
else
    U=1;
    for j=1:t,U=kron(U,V);end
end
end
