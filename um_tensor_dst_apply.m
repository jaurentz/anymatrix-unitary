function Y=um_tensor_dst_apply(X,grid)
%UM_TENSOR_DST_APPLY Separable orthonormal DST-I on 2D/3D tensor grids.
%   Each dimension uses shared O(n log n) um_dst1_apply. GRID(1)
%   is the fastest physical index in MATLAB column-major storage.
%   The input/output have prod(GRID) rows and any number of RHS columns.
d=numel(grid);N=prod(grid);
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_tensor_dst_apply:Dimension', ...
        'Input must have prod(GRID) rows.');
end
if issparse(X),X=full(X);end
nrhs=size(X,2);
Y=reshape(X,[grid,nrhs]);
for j=1:d
    other=setdiff(1:d,j,'stable');
    permutation=[j,other,d+1];
    W=permute(Y,permutation);
    W=um_dst1_apply(reshape(W,grid(j),[]));
    W=reshape(W,[grid(permutation(1:d)),nrhs]);
    Y=ipermute(W,permutation);
end
Y=reshape(Y,N,nrhs);
end
