function Y=um_tensor_permute_apply(X,flag,d,pi)
%UM_TENSOR_PERMUTE_APPLY Reshape/permute tensor factors, no matrix assembly.
%   PI maps (i_1,...,i_t) to (i_pi(1),...,i_pi(t)) in kron order.
%   MATLAB uses column-major reshapes, hence reversed axis indexing.
t=numel(pi);
n=d^t;
if strcmpi(flag,'size'),Y=[n n];return;end
if any(strcmpi(flag,{'transp','ctransp'}))
    invpi=zeros(size(pi));
    invpi(pi)=1:t;
    pi=invpi;
elseif ~strcmpi(flag,'notransp')
    error('um_tensor_permute_apply:Flag','Invalid operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_tensor_permute_apply:Dimension','Input must have d^t rows.');
end
if t==1,Y=X;return;end
axes=t+1-pi(t:-1:1);
Y=full(X);
for j=1:size(Y,2)
    B=reshape(Y(:,j),d*ones(1,t));
    Y(:,j)=reshape(permute(B,axes),n,1);
end
end
