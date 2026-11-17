function [U,properties] = haar_tensor_power(d,t,seed)
%HAAR_TENSOR_POWER The same Haar U(d) repeated on t tensor factors.
%   U=kron(V,...,V), t copies, V=haar_unitary(d,seed).
%   Eigenvalues are all ordered t-fold products of eigenvalues of V.
%   tr(U^k)/d^t = (tr(V^k)/d)^t, giving exact reference moments.
%   This is the FULL tensor space, not a symmetric or exterior power.
%   Ref: L. Marcinkowska et al., arXiv:2605.31168 (2026).
%   https://arxiv.org/abs/2605.31168
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(t),t=2;end
validateattributes(t,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3,seed=[];end
if d^t>4096
    error('haar_tensor_power:Order','Explicit order >4096; reduce d or t.');
end
V=haar_unitary(d,seed);
U=1;
for j=1:t,U=kron(U,V);end
end
