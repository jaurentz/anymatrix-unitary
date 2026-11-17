function [U,properties] = haar_block_tensor_power(d,t1,t2,seed,output)
%HAAR_BLOCK_TENSOR_POWER Two independent Haar base factors on tensor blocks.
%   U=V1^(tensor t1) kron V2^(tensor t2), with independent V1,V2 in U(d).
%   OUTPUT='operator' avoids explicit Kronecker matrices, supports multi-RHS
%   forward and adjoint actions and reports size with U([],'size').
%   OUTPUT='spectral' requires a SEED, returning exact finite sample
%   eigenvalues/moments computed only from the two small base matrices.
%   Normalized moments: (tr(V1^k)/d)^t1*(tr(V2^k)/d)^t2.
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
if nargin<5||isempty(output),output='matrix';end
n=d^(t1+t2);
if ~isfinite(n)||n>flintmax
    error('haar_block_tensor_power:Order','Tensor order exceeds exact index range.');
end
if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('haar_block_tensor_power',d,t1,t2,seed);
    return
end
implicit=um_is_operator(output);
if ~implicit && n>4096
    error('haar_block_tensor_power:Order','Explicit order >4096; use operator mode.');
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
if implicit
    factors=[repmat({V1},1,t1),repmat({V2},1,t2)];
    U=@(X,flag)um_tensor_apply(X,flag,factors);
else
    U=1;
    for j=1:t1,U=kron(U,V1);end
    for j=1:t2,U=kron(U,V2);end
end
end
