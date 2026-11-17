function [U,properties] = tensor_permutation(d,pi)
%TENSOR_PERMUTATION Permute factors of a d^t-dimensional tensor product.
%   PI is a permutation of 1:t. Basis column (i1,...,it) maps to
%   (i_PI(1),...,i_PI(t)), in MATLAB kron/big-endian ordering.
%   U is a REAL permutation matrix, hence orthogonal and unitary.
%   For any k, tr(U^k)=d^(number of cycles of PI^k), an exact formula.
%   Special case PI=[2 1] is the swap matrix with known +/-1 counts.
%   Schur-Weyl duality and quantum channel twirling motivate this family.
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
U=sparse(row+1,cols+1,1,n,n);
end
