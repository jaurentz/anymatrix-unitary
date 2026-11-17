function row=um_tensor_permutation_map(d,pi)
%UM_TENSOR_PERMUTATION_MAP Exact 1-based row index for each basis column.
%   Row = permuted tensor-index tuple in Kronecker (big-endian) order.
t=numel(pi);
n=d^t;
if ~isfinite(n)||n>2e6
    error('tensor_permutation:Order','Explicit order exceeds 2e6.');
end
cols=(0:n-1)';
digits=zeros(n,t);
for j=1:t
    digits(:,j)=mod(floor(cols/d^(t-j)),d);
end
row=zeros(n,1);
for j=1:t
    row=row+digits(:,pi(j))*d^(t-j);
end
row=row+1;
end
