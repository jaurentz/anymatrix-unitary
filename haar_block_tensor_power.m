function [U,properties] = haar_block_tensor_power(d,t1,t2,seed)
%HAAR_BLOCK_TENSOR_POWER Independent Haar factors reused within two blocks.
%   U=V1^(tensor t1) kron V2^(tensor t2), where V1,V2 are independent
%   Haar U(d) samples. T1,T2>=1. Compare with HAAR_TENSOR_POWER,
%   which reuses the SAME unitary in all tensor factors.
%   Exact moment: [tr(V1^k)/d]^t1*[tr(V2^k)/d]^t2.
%   Ref: L. Marcinkowska et al., arXiv:2605.31168 (2026).
%   https://arxiv.org/abs/2605.31168
properties={'square','unitary','normal','complex','scalable','random'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(t1),t1=1;end
if nargin<3||isempty(t2),t2=1;end
validateattributes(t1,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(t2,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<4,seed=[];end
if d^(t1+t2)>4096
    error('haar_block_tensor_power:Order','Explicit order >4096; reduce factors.');
end
if ~isempty(seed)
    validateattributes(seed,{'numeric'}, ...
        {'scalar','integer','nonnegative','finite'});
    state=rng;
    cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
    rng(seed,'twister');
end
V1=haar_unitary(d);
V2=haar_unitary(d);
U=1;
for j=1:t1,U=kron(U,V1);end
for j=1:t2,U=kron(U,V2);end
end
