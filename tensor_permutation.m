function [U,properties] = tensor_permutation(d,pi,output)
%TENSOR_PERMUTATION Orthogonal representation of a tensor-factor permutation.
%   PI is a permutation of 1:t. U maps basis (i1,...,it) to
%   (i_PI(1),...,i_PI(t)) in Kronecker (big-endian) order.
%   OUTPUT='operator' permutes reshaped tensors with no sparse matrix.
%   OUTPUT='spectral' provides exact trace moments from permutation cycles,
%   with eigenphases obtainable from basis orbits for moderate dimensions.
%   Trace(U^k)=d^cycles(PI^k); swap PI=[2 1] has only +/-1 eigenvalues.
%   Ref: https://doi.org/10.1088/1751-8121/ae965b
%   https://arxiv.org/abs/2605.31168
properties={'square','unitary','normal','real','orthogonal','sparse', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(pi),pi=[2 1];end
validateattributes(pi,{'numeric'},{'vector','integer','positive','finite'});
pi=pi(:)';
t=numel(pi);
if ~isequal(sort(pi),1:t)
    error('tensor_permutation:Permutation','PI must permute 1:t.');
end
if nargin<3||isempty(output),output='matrix';end
n=d^t;
if ~isfinite(n)||n>flintmax
    error('tensor_permutation:Order','d^t exceeds exact index range.');
end
if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('tensor_permutation',d,pi);
    return
end
if um_is_operator(output)
    U=@(X,flag)um_tensor_permute_apply(X,flag,d,pi);
else
    row=um_tensor_permutation_map(d,pi);
    U=sparse(row,1:n,1,n,n);
end
end
